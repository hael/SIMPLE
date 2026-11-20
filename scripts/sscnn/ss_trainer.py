#!/usr/bin/env python3
"""Secondary-structure selection trainer for 2D class averages.

You label class averages in a browser (helices visible: yes / no / unsure); a logistic model learns from your yes/no
labels and predicts the others; you correct and retrain; the selection is exported as a stack of the chosen class
averages for `simple_exec prg=map_cavgs_selection`.

The model sees six measurements per class, chosen for helices (rods about 10 A wide, about 10 A apart), each
computed against the noise of the even/odd half-averages, so that a feature counts only when both halves carry it:
  frc_6_12      mean FRC of the even and odd halves over 12-6 A (signal in the helix band at all)
  ridge_snr     rod-like ridge strength (Hessian at a 3.5 A scale) of the average over that of the half difference
  band_ratio    signal power at 12-8 A over that at 30-15 A (helix spacing against the overall shape)
  band_aniso    angular anisotropy of the 12-8 A signal power (parallel helices in side views)
  log_pop       log10 of the class population
  neg_log_res   minus log of the class resolution (if a cls2D text is given)
With 50-100 labels this is what can be learned; a CNN needs the simulated training data of a later step.

Usage (numpy required):
  ss_trainer.py --stk cavgs_iterNNN.mrc --mskdiam A --workdir DIR [--cls cls2D.txt] [--port 8765]
then open http://127.0.0.1:8765. The even/odd halves are read from <stk>_even.mrc and <stk>_odd.mrc.
Labels persist in DIR/labels.json; the export writes DIR/selected.mrc and DIR/selected_classes.txt.
"""
import argparse
import base64
import json
import os
import re
import struct
import sys
import threading
import zlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

FEATURES = ['frc_6_12', 'ridge_snr', 'band_ratio', 'band_aniso', 'log_pop', 'neg_log_res']
RIDGE_SIGMA_A = 3.5
L2 = 1.0


# ------------------------------------------------------------------ data
def read_mrc(path):
    with open(path, 'rb') as f:
        h = f.read(1024)
        nx, ny, nz, mode = struct.unpack('<4i', h[:16])
        xlen = struct.unpack('<f', h[40:44])[0]
        nsym = struct.unpack('<i', h[92:96])[0]
        if mode != 2:
            raise SystemExit(f'{path}: MRC mode {mode}, expected 2 (32-bit reals)')
        f.seek(1024 + nsym)
        a = np.fromfile(f, dtype='<f4', count=nx * ny * nz).reshape(nz, ny, nx)
    return a, xlen / nx, h, nsym


def write_selected(src, inds, out):
    """The images inds (0-based) of src copied byte for byte into a new stack, the header's count updated."""
    _, smpd, h, nsym = read_mrc(src)
    nx, ny = struct.unpack('<2i', h[:8])
    hb = bytearray(h)
    struct.pack_into('<i', hb, 8, len(inds))
    struct.pack_into('<f', hb, 48, len(inds) * smpd)      # the cell's z length
    struct.pack_into('<i', hb, 92, 0)                     # no extended header
    nb = nx * ny * 4
    with open(src, 'rb') as f, open(out, 'wb') as g:
        g.write(bytes(hb))
        for i in inds:
            f.seek(1024 + nsym + i * nb)
            g.write(f.read(nb))


def read_cls(path, n):
    pat = {k: re.compile(r'(?:^| )' + k + r'=(\S+)') for k in ('pop', 'res')}
    pop, res = np.full(n, np.nan), np.full(n, np.nan)
    if path:
        for i, line in enumerate(open(path)):
            if i >= n:
                break
            for k, arr in (('pop', pop), ('res', res)):
                m = pat[k].search(line)
                if m:
                    arr[i] = float(m.group(1))
    return pop, res


# ------------------------------------------------------------------ measurements
class Measurer:
    def __init__(self, box, smpd, mskdiam):
        self.box, self.smpd = box, smpd
        c = box // 2
        y, x = np.mgrid[0:box, 0:box] - c
        self.r = np.hypot(x, y) * smpd
        rad = mskdiam / 2.0
        self.rad = rad
        self.mask = np.clip((rad - self.r) / 10.0 + 0.5, 0.0, 1.0)
        self.inner = self.r <= rad
        fx = np.fft.fftfreq(box, d=smpd)
        FX, FY = np.meshgrid(fx, fx)
        self.F = np.hypot(FX, FY)
        self.ang = np.arctan2(FY, FX) % np.pi
        g = np.exp(-2 * np.pi ** 2 * RIDGE_SIGMA_A ** 2 * self.F ** 2) * RIDGE_SIGMA_A ** 2
        self.hxx = -(2 * np.pi * FX) ** 2 * g
        self.hyy = -(2 * np.pi * FY) ** 2 * g
        self.hxy = -(2 * np.pi) ** 2 * FX * FY * g
        self.shell = np.rint(self.F * box * smpd).astype(int)
        self.band_hx = (self.F >= 1 / 12.0) & (self.F <= 1 / 6.0)
        self.band_h = (self.F >= 1 / 12.0) & (self.F <= 1 / 8.0)
        self.band_l = (self.F >= 1 / 30.0) & (self.F <= 1 / 15.0)

    def prep(self, im):
        out = self.r > self.rad
        bg = im[out].mean() if out.any() else im.mean()
        return (im - bg) * self.mask

    def ridge(self, ft):
        hxx = np.real(np.fft.ifft2(ft * self.hxx))
        hyy = np.real(np.fft.ifft2(ft * self.hyy))
        hxy = np.real(np.fft.ifft2(ft * self.hxy))
        tr, d = 0.5 * (hxx + hyy), np.sqrt(0.25 * (hxx - hyy) ** 2 + hxy ** 2)
        l1, l2 = tr + d, tr - d
        big = np.maximum(np.abs(l1), np.abs(l2))
        small = np.minimum(np.abs(l1), np.abs(l2))
        return (big - small)[self.inner].mean()

    def measure(self, avg, ev, od):
        a, e, o = self.prep(avg), self.prep(ev), self.prep(od)
        A, E, O = np.fft.fft2(a), np.fft.fft2(e), np.fft.fft2(o)
        D = 0.5 * (E - O)                                   # the noise of the average, without its signal
        # FRC of the halves over the helix band, shell by shell
        shells = np.unique(self.shell[self.band_hx])
        frc = []
        for s in shells:
            m = self.shell == s
            num = np.real((E[m] * np.conj(O[m])).sum())
            den = np.sqrt((np.abs(E[m]) ** 2).sum() * (np.abs(O[m]) ** 2).sum())
            if den > 0:
                frc.append(num / den)
        frc_band = float(np.mean(frc)) if frc else 0.0
        rn = self.ridge(D)
        ridge_snr = self.ridge(A) / rn if rn > 0 else 0.0
        pa, pd = np.abs(A) ** 2, np.abs(D) ** 2
        sig_h = max(pa[self.band_h].sum() - pd[self.band_h].sum(), 0.0)
        sig_l = max(pa[self.band_l].sum() - pd[self.band_l].sum(), 1e-12)
        band_ratio = sig_h / sig_l
        # anisotropy: the order-2 harmonic of the 12-8 A signal power over its mean
        w = np.clip(pa - pd, 0, None)[self.band_h]
        th = self.ang[self.band_h]
        tot = w.sum()
        aniso = float(np.abs((w * np.exp(2j * th)).sum()) / tot) if tot > 0 else 0.0
        return [frc_band, ridge_snr, band_ratio, aniso]


# ------------------------------------------------------------------ model
def fit_logistic(X, y, l2=L2, iters=50):
    """Logistic regression by Newton steps with an L2 penalty (not on the intercept); returns weights with intercept."""
    Xb = np.hstack([np.ones((len(X), 1)), X])
    w = np.zeros(Xb.shape[1])
    P = np.eye(Xb.shape[1]) * l2
    P[0, 0] = 0.0
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-Xb @ w))
        g = Xb.T @ (p - y) + P @ w
        H = Xb.T @ (Xb * (p * (1 - p))[:, None]) + P + 1e-9 * np.eye(len(w))
        step = np.linalg.solve(H, g)
        w -= step
        if np.abs(step).max() < 1e-8:
            break
    return w


def predict(w, X):
    return 1.0 / (1.0 + np.exp(-(np.hstack([np.ones((len(X), 1)), X]) @ w)))


def auc(pos, neg):
    if len(pos) == 0 or len(neg) == 0:
        return float('nan')
    pos, neg = np.asarray(pos), np.asarray(neg)
    return float(((pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()) / (len(pos) * len(neg)))


# ------------------------------------------------------------------ images for the page
def png(img, size):
    n = img.shape[0]
    idx = (np.arange(size) * n / size).astype(int)
    im = img[np.ix_(idx, idx)]
    m, s = im.mean(), im.std() or 1.0
    g = np.clip((im - (m - 3 * s)) / (6 * s) * 255, 0, 255).astype(np.uint8)
    raw = b''.join(b'\x00' + g[j].tobytes() for j in range(size))

    def chunk(tag, data):
        c = tag + data
        return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)
    return (b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', size, size, 8, 0, 0, 0, 0))
            + chunk(b'IDAT', zlib.compress(raw, 6)) + chunk(b'IEND', b''))


# ------------------------------------------------------------------ the trainer
class Trainer:
    def __init__(self, a):
        self.stk, self.workdir = a.stk, a.workdir
        os.makedirs(a.workdir, exist_ok=True)
        self.avg, self.smpd, _, _ = read_mrc(a.stk)
        base = a.stk[:-4] if a.stk.endswith('.mrc') else a.stk
        ev, od = base + '_even.mrc', base + '_odd.mrc'
        if not (os.path.exists(ev) and os.path.exists(od)):
            raise SystemExit(f'the half-averages {ev} and {od} are needed')
        self.even, _, _, _ = read_mrc(ev)
        self.odd, _, _, _ = read_mrc(od)
        n, box = self.avg.shape[0], self.avg.shape[1]
        self.n = n
        pop, res = read_cls(a.cls, n)
        self.pop, self.res = pop, res
        self.empty = np.array([not np.isfinite(self.avg[i]).all() or self.avg[i].std() == 0 or
                               (np.isfinite(pop[i]) and pop[i] <= 0) for i in range(n)])
        cache = os.path.join(a.workdir, 'features.npy')
        if os.path.exists(cache) and np.load(cache).shape == (n, len(FEATURES)):
            self.X = np.load(cache)
        else:
            m = Measurer(box, self.smpd, a.mskdiam)
            X = np.full((n, len(FEATURES)), np.nan)
            for i in range(n):
                if self.empty[i]:
                    continue
                X[i, :4] = m.measure(self.avg[i], self.even[i], self.odd[i])
                X[i, 4] = np.log10(pop[i]) if np.isfinite(pop[i]) and pop[i] > 0 else np.nan
                X[i, 5] = -np.log(res[i]) if np.isfinite(res[i]) and res[i] > 0 else np.nan
            np.save(cache, X)
            self.X = X
        self.cols = [k for k in range(len(FEATURES)) if np.isfinite(self.X[~self.empty, k]).all()]
        self.labels_path = os.path.join(a.workdir, 'labels.json')
        self.labels = json.load(open(self.labels_path)) if os.path.exists(self.labels_path) else {}
        self.p = np.full(n, np.nan)
        self.metrics = 'not trained yet'
        self.lock = threading.Lock()
        self.train()

    def save_labels(self):
        tmp = self.labels_path + '.tmp'
        json.dump(self.labels, open(tmp, 'w'), indent=0, sort_keys=True)
        os.replace(tmp, self.labels_path)

    def train(self):
        live = np.where(~self.empty)[0]
        Z = self.X[:, self.cols]
        mu, sd = Z[live].mean(0), Z[live].std(0)
        sd[sd == 0] = 1.0
        Z = (Z - mu) / sd
        lab = [(int(k) - 1, v) for k, v in self.labels.items() if v in ('yes', 'no')]
        lab = [(i, v) for i, v in lab if 0 <= i < self.n and not self.empty[i]]
        ny = sum(v == 'yes' for _, v in lab)
        nn = len(lab) - ny
        if ny < 3 or nn < 3:
            self.p[:] = np.nan
            self.metrics = f'needs at least 3 "yes" and 3 "no" labels (has {ny} and {nn})'
            return
        idx = np.array([i for i, _ in lab])
        y = np.array([1.0 if v == 'yes' else 0.0 for _, v in lab])
        w = fit_logistic(Z[idx], y)
        self.p[:] = np.nan
        self.p[live] = predict(w, Z[live])
        loo = np.array([predict(fit_logistic(np.delete(Z[idx], q, 0), np.delete(y, q)), Z[idx[q:q + 1]])[0]
                        for q in range(len(idx))])
        acc = ((loo >= 0.5) == (y == 1)).mean()
        coef = ', '.join(f'{FEATURES[c]} {w[q + 1]:+.2f}' for q, c in enumerate(self.cols))
        self.metrics = (f'{ny} yes, {nn} no labels; leave-one-out accuracy {100 * acc:.0f}%, '
                        f'AUC {auc(loo[y == 1], loo[y == 0]):.2f}; weights: {coef}')

    def state(self):
        rows = []
        for i in range(self.n):
            rows.append({'cls': i + 1, 'empty': bool(self.empty[i]),
                         'pop': None if not np.isfinite(self.pop[i]) else int(self.pop[i]),
                         'res': None if not np.isfinite(self.res[i]) else round(float(self.res[i]), 1),
                         'label': self.labels.get(str(i + 1), ''),
                         'p': None if not np.isfinite(self.p[i]) else round(float(self.p[i]), 3),
                         'feat': {FEATURES[k]: (None if not np.isfinite(self.X[i, k]) else round(float(self.X[i, k]), 3))
                                  for k in range(len(FEATURES))}})
        return {'rows': rows, 'metrics': self.metrics, 'stk': self.stk, 'smpd': round(self.smpd, 4)}

    def export(self, thr):
        sel = []
        for i in range(self.n):
            if self.empty[i]:
                continue
            lab = self.labels.get(str(i + 1), '')
            if lab == 'yes' or (lab != 'no' and np.isfinite(self.p[i]) and self.p[i] >= thr):
                sel.append(i)
        if not sel:
            return 'nothing selected'
        out = os.path.join(self.workdir, 'selected.mrc')
        write_selected(self.stk, sel, out)
        with open(os.path.join(self.workdir, 'selected_classes.txt'), 'w') as f:
            f.write('# classes with secondary structure: labelled yes, or unlabelled with p >= %.2f\n' % thr)
            f.write('\n'.join(str(i + 1) for i in sel) + '\n')
        return (f'{len(sel)} classes written to {out} (list in selected_classes.txt). To apply to a project:\n'
                f'simple_exec prg=map_cavgs_selection projfile=<project>.simple stk={self.stk} stk2={out}')


PAGE = r"""<!doctype html><meta charset="utf-8"><title>Secondary-structure selection</title>
<style>
body{background:#fcfcfb;color:#0b0b0b;font:13px system-ui,sans-serif;margin:16px}
#bar{position:sticky;top:0;background:#fcfcfb;padding:6px 0 8px;border-bottom:1px solid #e4e3df;z-index:2}
button,select{font:inherit;padding:3px 8px}
#grid{display:flex;flex-wrap:wrap;gap:8px;margin-top:10px}
.t{width:136px;cursor:pointer}.t img{width:128px;height:128px;display:block;border:4px solid #e4e3df;border-radius:4px}
.t.yes img{border-color:#2a78d6}.t.no img{border-color:#eb6834}.t.unsure img{border-color:#eda100}
.t.empty{opacity:.35}.t.focus img{outline:2px solid #0b0b0b;outline-offset:1px}
.cap{font-size:11px;color:#52514e;line-height:1.3}.pbar{height:4px;background:#e4e3df;border-radius:2px;margin-top:2px}
.pbar i{display:block;height:4px;background:#2a78d6;border-radius:2px}
#big{position:fixed;inset:0;background:rgba(252,252,251,.97);display:none;z-index:5;padding:20px}
#big img{width:384px;height:384px;border:1px solid #e4e3df;margin-right:8px}
#msg{white-space:pre-wrap;color:#52514e;margin-top:4px}
.key{display:inline-block;border:1px solid #c3c2b7;border-radius:3px;padding:0 4px;font-size:11px}
</style>
<div id="bar">
 <b>Helices visible?</b> Hover or click a class, then <span class="key">y</span> yes,
 <span class="key">n</span> no, <span class="key">u</span> unsure, <span class="key">0</span> clear;
 <span class="key">Enter</span> enlarges (average, even, odd), <span class="key">←</span><span class="key">→</span> move, <span class="key">Esc</span> closes.
 &nbsp; Sort <select id="sort"><option value="cls">class</option><option value="unc">most uncertain first</option>
 <option value="p">predicted yes first</option><option value="res">resolution</option><option value="pop">population</option></select>
 &nbsp;<button id="train">Train</button>
 &nbsp;Export threshold <input id="thr" type="number" value="0.5" min="0" max="1" step="0.05" style="width:4em">
 <button id="export">Export selection</button>
 &nbsp;<span style="color:#2a78d6">■</span> yes <span style="color:#eb6834">■</span> no <span style="color:#eda100">■</span> unsure; bar: predicted p(yes)
 <div id="msg"></div>
</div>
<div id="grid"></div>
<div id="big"><div id="bigcap" style="margin-bottom:8px"></div><img id="b0"><img id="b1"><img id="b2">
<div class="cap">average · even half · odd half — keys work here too</div></div>
<script>
let S=/*INIT*/null, order=[], focus=-1;
const $=q=>document.querySelector(q);
async function load(){S=await (await fetch('/state')).json(); $('#msg').textContent=S.metrics; render();}
function sorted(){const r=S.rows.slice(), k=$('#sort').value;
 const key={cls:x=>x.cls,unc:x=>x.p==null?9:Math.abs(x.p-.5),p:x=>x.p==null?9:-x.p,res:x=>x.res??1e9,pop:x=>-(x.pop??0)}[k];
 r.sort((a,b)=>(a.empty-b.empty)||(key(a)-key(b))); return r;}
function cap(x){return `class ${x.cls} · n=${x.pop??'?'} · ${x.res??'?'} Å`+(x.p==null?'':` · p=${x.p.toFixed(2)}`);}
function render(){order=sorted(); const g=$('#grid'); g.innerHTML='';
 order.forEach((x,q)=>{const d=document.createElement('div'); d.className='t '+(x.label||'')+(x.empty?' empty':'')+(q==focus?' focus':'');
  const f=Object.entries(x.feat).map(([k,v])=>`${k} ${v}`).join('\n');
  d.innerHTML=`<img src="/img/${x.cls}/avg/128" title="${f}"><div class="cap">${cap(x)}</div>`+
   (x.p==null?'':`<div class="pbar"><i style="width:${(x.p*100).toFixed(0)}%"></i></div>`);
  d.onmouseenter=()=>setFocus(q,false); d.onclick=()=>{setFocus(q,false); openBig();}; g.appendChild(d);});}
function setFocus(q,scroll){const t=document.querySelectorAll('.t'); if(focus>=0&&t[focus])t[focus].classList.remove('focus');
 focus=q; if(t[q]){t[q].classList.add('focus'); if(scroll)t[q].scrollIntoView({block:'nearest'});} if($('#big').style.display=='block')fillBig();}
function fillBig(){const x=order[focus]; $('#bigcap').textContent=cap(x)+(x.label?` · labelled ${x.label}`:'');
 ['avg','even','odd'].forEach((k,q)=>$('#b'+q).src=`/img/${x.cls}/${k}/384`);}
function openBig(){if(focus<0)return; $('#big').style.display='block'; fillBig();}
async function setLabel(v){const x=order[focus]; if(!x||x.empty)return;
 await fetch('/label',{method:'POST',body:JSON.stringify({cls:x.cls,label:v})}); x.label=v;
 const t=document.querySelectorAll('.t')[focus]; t.className='t '+v+' focus'; if($('#big').style.display=='block')fillBig();
 setFocus(Math.min(focus+1,order.length-1),true);}
document.addEventListener('keydown',e=>{if(e.target.tagName=='INPUT')return;
 if(e.key=='y')setLabel('yes'); else if(e.key=='n')setLabel('no'); else if(e.key=='u')setLabel('unsure'); else if(e.key=='0')setLabel('');
 else if(e.key=='ArrowRight')setFocus(Math.min(focus+1,order.length-1),true); else if(e.key=='ArrowLeft')setFocus(Math.max(focus-1,0),true);
 else if(e.key=='Enter')openBig(); else if(e.key=='Escape')$('#big').style.display='none';});
$('#sort').onchange=()=>{focus=-1; render();};
$('#train').onclick=async()=>{$('#msg').textContent='training…'; await fetch('/train',{method:'POST'}); focus=-1; load();};
$('#export').onclick=async()=>{const r=await (await fetch('/export',{method:'POST',body:JSON.stringify({thr:+$('#thr').value})})).json(); $('#msg').textContent=r.message;};
if(S){$('#msg').textContent=S.metrics; render();} else load();
</script>"""


def make_handler(tr):
    class H(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def send(self, body, ctype, code=200):
            self.send_response(code)
            self.send_header('Content-Type', ctype)
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path == '/':
                with tr.lock:
                    init = json.dumps(tr.state()).replace('</', '<\\/')
                return self.send(PAGE.replace('/*INIT*/null', init).encode(), 'text/html; charset=utf-8')
            if self.path == '/state':
                with tr.lock:
                    return self.send(json.dumps(tr.state()).encode(), 'application/json')
            m = re.fullmatch(r'/img/(\d+)/(avg|even|odd)/(\d+)', self.path)
            if m:
                i, kind, size = int(m.group(1)) - 1, m.group(2), min(int(m.group(3)), 512)
                if 0 <= i < tr.n:
                    src = {'avg': tr.avg, 'even': tr.even, 'odd': tr.odd}[kind]
                    return self.send(png(src[i], size), 'image/png')
            self.send(b'not found', 'text/plain', 404)

        def do_POST(self):
            body = self.rfile.read(int(self.headers.get('Content-Length', 0)) or 0)
            data = json.loads(body) if body else {}
            with tr.lock:
                if self.path == '/label':
                    k, v = str(int(data['cls'])), data.get('label', '')
                    if v in ('yes', 'no', 'unsure'):
                        tr.labels[k] = v
                    else:
                        tr.labels.pop(k, None)
                    tr.save_labels()
                    return self.send(b'{}', 'application/json')
                if self.path == '/train':
                    tr.train()
                    return self.send(json.dumps({'metrics': tr.metrics}).encode(), 'application/json')
                if self.path == '/export':
                    msg = tr.export(float(data.get('thr', 0.5)))
                    return self.send(json.dumps({'message': msg}).encode(), 'application/json')
            self.send(b'not found', 'text/plain', 404)
    return H


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stk', required=True, help='class averages (MRC); <stk>_even.mrc and <stk>_odd.mrc beside it')
    ap.add_argument('--mskdiam', type=float, required=True, help='mask diameter (A) the measurements are taken in')
    ap.add_argument('--workdir', required=True, help='where labels, features and the export go')
    ap.add_argument('--cls', default='', help='cls2D text (simple_exec prg=print_project_field oritype=cls2D)')
    ap.add_argument('--port', type=int, default=8765)
    a = ap.parse_args()
    tr = Trainer(a)
    print(f'{tr.n} class averages, {int((~tr.empty).sum())} with content; {tr.metrics}')
    print(f'open http://127.0.0.1:{a.port}  (Ctrl-C stops; labels are saved as you go)')
    sys.stdout.flush()
    ThreadingHTTPServer(('127.0.0.1', a.port), make_handler(tr)).serve_forever()


if __name__ == '__main__':
    main()
