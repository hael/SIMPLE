#!/usr/bin/env python3
"""A labelling page for real class averages, to train the class-average CNN on them. Each class gets one of
  secondary_structure  protein with helices or sheets visible
  likely_protein       protein, no secondary structure visible
  protein              protein, unsure whether secondary structure shows
  junk                 ice, carbon, aggregates, contamination, noise
  unsure               cannot tell protein from junk
on a self-contained HTML page that needs no server: make it where the data is, label in any browser, export a TSV.

  python3 label_cavgs.py PATH [PATH ...] --out labels.html [--dataset NAME] [--blind] [--simple SIMPLE_BUILD]

PATH, any mix of:
  - a project with class averages (a .simple file, e.g. solve2D's): its classes in class order, without CNN calls;
  - a cavg_cnn_predictions.txt (from simple_exec prg=cavg_cnn_predict): the classes of the project copy beside it,
    with the CNN's calls;
  - a folder: every cavg_cnn_predictions.txt under it (a cavg_cnn_predict execution folder, a mic_training_generator.py output).
Classes without particles are left out, and class averages already on the page are not shown twice. --dataset names
the data (default: the project's name); the training's cross-validation will hold out whole datasets, so give runs
of the same sample the same name. --blind leaves the CNN's calls out (tiles in class order, no probabilities), so
they cannot sway the labels; otherwise the tiles of a table are grouped by the CNN's call with its probability.

On the page: pick a label (its button, or the keys s l p j u; 0 clears), then click tiles to give it; clicking a tile
that has it clears it. 'Label the rest' gives it to every unlabelled tile of a section. Tile size, brightness and
contrast are sliders. Labels are kept in the browser (localStorage) as you go; 'Export TSV' downloads them and 'Import
TSV' loads an export to continue or revise. The TSV has one row per class shown: dataset, run, stk, source (the
project or table), cls, pop, res, label (empty: not labelled), cnn_group (empty without the CNN's call), blind,
curator, time.
Standard library only; SIMPLE (--simple, else $SIMPLE_PATH, else the install or build this script sits in) reads
the projects.
"""
import argparse
import base64
import hashlib
import html
import importlib.util
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('mic_training_generator', os.path.join(HERE, '..', 'mic_training_generator.py'))
MS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(MS)                                  # read_mrc, png, records, find_simple and the colours

TABLE = 'cavg_cnn_predictions.txt'
CNN_GROUPS = ('secondary_structure', 'likely_protein', 'junk')
PROB = {'secondary_structure': 'p_ss', 'likely_protein': 'p_likely', 'junk': 'p_junk'}
MAX_PX = 256                                                  # images are embedded at their box, at most this


def natural(s):
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', s)]


def find_inputs(paths):
    """(kind, file, root) per input: a project ('.simple'), a prediction table, or every table under a folder, in
    natural order."""
    found = []
    for p in paths:
        p = os.path.abspath(p)
        if os.path.isfile(p):
            found.append(('project' if p.endswith('.simple') else 'table', p, os.path.dirname(p)))
        elif os.path.isdir(p):
            hits = [os.path.join(r, TABLE) for r, _, files in os.walk(p) if TABLE in files]
            if not hits:
                raise SystemExit(f'no {TABLE} under {p}: give a project as its .simple file')
            found += [('table', t, p) for t in sorted(hits, key=lambda t: natural(os.path.relpath(t, p)))]
        else:
            raise SystemExit(f'no file or folder {p}')
    return found


def run_name(kind, path, root):
    """A table's run: its folder, or the folder above a numbered cavg_cnn_predict folder, relative to the folder
    above root. A project's: its folder and the one above."""
    d = os.path.dirname(path)
    if kind == 'project':
        return os.path.join(os.path.basename(os.path.dirname(d)), os.path.basename(d))
    if re.match(r'\d+_cavg_cnn_predict$', os.path.basename(d)):
        d = os.path.dirname(d)
    r = os.path.relpath(d, os.path.dirname(root))
    return os.path.basename(d) if r == '.' or r.startswith('..') else r


def cavgs_of(exe, env, proj):
    """The project's class averages, chosen as SIMPLE's get_cavgs_stk does: the one out entry of imgkind cavg."""
    outs = [o for o in MS.records(exe, env, proj, 'out', os.path.dirname(proj)) if o.get('imgkind') == 'cavg']
    if len(outs) != 1:
        raise SystemExit(f'{proj}: {len(outs)} class-average entries in the project, expected one')
    stk = outs[0].get('stk', '')
    if not os.path.exists(stk) and 'stkpath' in outs[0]:
        stk = os.path.join(outs[0]['stkpath'], stk)
    if not os.path.exists(stk):
        raise SystemExit(f'{proj}: its class averages {stk} do not exist')
    return os.path.abspath(stk)


def res_of(r):
    """The class's resolution, None when absent or not finite (JSON has no NaN)."""
    try:
        v = float(r.get('res', ''))
    except ValueError:
        return None
    return round(v, 2) if math.isfinite(v) else None


def collect(a, exe, env):
    items, sections, seen = [], [], set()
    for kind, path, root in find_inputs(a.paths):
        if kind == 'table':                                   # the project copy cavg_cnn_predict made beside it
            folder = os.path.dirname(path)
            projs = [f for f in os.listdir(folder) if f.endswith('.simple')]
            if len(projs) != 1:
                raise SystemExit(f'{folder}: expected one project beside {TABLE}, found {len(projs)}')
            proj = os.path.join(folder, projs[0])
            with open(path) as f:
                rows = [dict(re.findall(r'(\w+)=(\S+)', l)) for l in f if '=' in l]
        else:                                                 # class, pop and res from cls2D; no CNN call
            proj = path
            rows = [{'class': str(int(float(r['class']))), 'pop': str(int(float(r.get('pop', 0)))),
                     'res': r.get('res', '')} for r in MS.records(exe, env, proj, 'cls2D', os.path.dirname(proj))]
            if not rows:
                raise SystemExit(f'{proj}: no classes (cls2D) in the project')
        stk = cavgs_of(exe, env, proj)
        imgs, box = MS.read_mrc(stk)
        if max(int(r['class']) for r in rows) > len(imgs):
            raise SystemExit(f'{stk}: {len(imgs)} images for {len(rows)} classes')
        size = min(box, MAX_PX)
        dataset = a.dataset or os.path.splitext(os.path.basename(proj))[0]
        run = run_name(kind, path, root)
        mine = []
        for r in rows:
            if int(r['pop']) == 0:
                continue
            key = f'{stk}#{r["class"]}'
            if key in seen:                                   # the same class averages from a second input
                continue
            seen.add(key)
            uri = 'data:image/png;base64,' + base64.b64encode(MS.png(imgs[int(r['class']) - 1], box, size)).decode()
            it = {'id': key, 'dataset': dataset, 'run': run, 'stk': stk, 'source': path, 'cls': int(r['class']),
                  'pop': int(r['pop']), 'res': res_of(r), 'img': uri}
            if kind == 'table' and not a.blind:
                it.update(cnn=r['group'], p=float(r[PROB[r['group']]]), pj=float(r['p_junk']), ps=float(r['p_ss']),
                          pl=float(r['p_likely']))
            mine.append(it)
        if not mine:
            continue
        if kind == 'project' or a.blind:
            groups = [{'title': '', 'ids': [it['id'] for it in mine]}]
        else:
            groups = []
            for g in CNN_GROUPS:
                sel = sorted((it for it in mine if it['cnn'] == g), key=lambda it: -it['p'])
                if sel:
                    groups.append({'title': f'CNN: {g.replace("_", " ")}', 'ids': [it['id'] for it in sel]})
        sections.append({'title': f'{dataset} · {run}: {len(mine)} classes', 'groups': groups})
        items += mine
        print(f'{run}: {len(mine)} classes from {stk}')
    return items, sections


STYLE = '''
:root{--ts:128px;--br:1;--ct:1}
body{background:%(s)s;color:%(t1)s;font:14px/1.45 system-ui,sans-serif;margin:0 24px 24px}
h1{font-size:20px;margin:16px 0 4px}h2{font-size:16px;margin:24px 0 4px}h3{font-size:14px;margin:14px 0 6px;color:%(t2)s}
p{color:%(t2)s;max-width:1100px}
.top{position:sticky;top:0;z-index:2;background:%(s)s;border-bottom:1px solid %(h)s}
.bar{padding:6px 0;display:flex;flex-wrap:wrap;gap:8px 16px;align-items:center;color:%(t2)s}
.bar label{display:flex;gap:6px;align-items:center}
button{font:inherit;padding:3px 9px;border:1px solid %(h)s;border-radius:4px;background:#fff;color:%(t1)s;cursor:pointer}
button.on{outline:2px solid %(t1)s}
button i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:-1px}
kbd{font:11px ui-monospace,monospace;color:%(t2)s;margin-left:4px}
.n{font-variant-numeric:tabular-nums;color:%(t2)s;margin-left:4px}
.g{display:flex;flex-wrap:wrap;gap:8px}
.t{width:var(--ts);position:relative;cursor:pointer}
.t img{display:block;width:var(--ts);height:var(--ts);border:3px solid %(h)s;border-radius:4px;box-sizing:border-box;
 filter:brightness(var(--br)) contrast(var(--ct))}
.t .b{position:absolute;top:5px;left:5px;min-width:16px;padding:0 3px;border-radius:3px;color:#fff;font:600 11px/16px
 system-ui,sans-serif;text-align:center}
.cap{font-size:11px;line-height:1.25;margin-top:2px;color:%(t2)s}
#msg{color:%(t1)s}
'''

SCRIPT = r'''
const D = JSON.parse(document.getElementById('data').textContent);
const L = {
  secondary_structure: {k: 's', s: 'S', c: '#2a6fdb', n: 'secondary structure'},
  likely_protein:      {k: 'l', s: 'L', c: '#1f9e74', n: 'likely protein'},
  protein:             {k: 'p', s: 'P', c: '#8a5cc2', n: 'protein, ss unsure'},
  junk:                {k: 'j', s: 'J', c: '#d4572a', n: 'junk'},
  unsure:              {k: 'u', s: 'U', c: '#6f6f6f', n: "can't tell"}};
const KEY = 'cavg_labels:' + D.page, CKEY = 'cavg_labels:curator';
const $ = id => document.getElementById(id);
let st = {}, store = true, brush = 'secondary_structure', dirty = false;
try { st = JSON.parse(localStorage.getItem(KEY) || '{}'); } catch (e) { store = false; }
const items = {}, tiles = {};
D.items.forEach(it => { items[it.id] = it; });

function save() {
  if (!store) return;
  try { localStorage.setItem(KEY, JSON.stringify(st)); } catch (e) { store = false; note(); }
}
function msg(s) { $('msg').textContent = s; }
function note() { if (!store) msg('labels are not kept by this browser: export before closing the page'); }
function fmt(v, d) { return v == null ? '?' : Number(v).toFixed(d); }

function paint(id) {
  const t = tiles[id], s = st[id], l = s ? s.l : '';
  t.img.style.borderColor = l ? L[l].c : '';
  t.b.textContent = l ? L[l].s : '';
  t.b.style.background = l ? L[l].c : 'transparent';
  t.b.title = l ? L[l].n : '';
}
function counts() {
  const n = {}; let tot = 0;
  Object.keys(L).forEach(l => { n[l] = 0; });
  D.items.forEach(it => { const s = st[it.id]; if (s && s.l in L) { n[s.l]++; tot++; } });
  Object.keys(L).forEach(l => { $('n_' + l).textContent = n[l]; });
  $('prog').textContent = `labelled ${tot} of ${D.items.length}`;
}
function set(id, l) {
  if (l) st[id] = {l: l, t: new Date().toISOString()}; else delete st[id];
  dirty = true; save(); paint(id); counts();
}
function pick(l) {
  brush = l;
  document.querySelectorAll('button[data-l]').forEach(b => b.classList.toggle('on', b.dataset.l === l));
}

function tile(it) {
  const d = document.createElement('div'); d.className = 't';
  const img = document.createElement('img'); img.src = it.img; img.loading = 'lazy';
  img.title = `class ${it.cls}, ${it.pop} particles, ${fmt(it.res, 1)} Å` +
    (it.cnn ? `; CNN junk ${fmt(it.pj, 3)}, ss ${fmt(it.ps, 3)}, likely ${fmt(it.pl, 3)}` : '');
  const b = document.createElement('span'); b.className = 'b';
  const c = document.createElement('div'); c.className = 'cap';
  c.textContent = `c${it.cls} · ${it.pop} ptcls · ${fmt(it.res, 1)} Å`;
  if (it.cnn) { const c2 = document.createElement('div'); c2.textContent = `CNN p ${fmt(it.p, 2)}`; c.appendChild(c2); }
  d.append(img, b, c);
  d.onclick = () => { const s = st[it.id]; set(it.id, brush === '' || (s && s.l === brush) ? '' : brush); };
  tiles[it.id] = {img: img, b: b};
  return d;
}

function build() {
  const bar = $('brushes');
  Object.keys(L).forEach(l => {
    const b = document.createElement('button'); b.dataset.l = l;
    b.innerHTML = `<i style="background:${L[l].c}"></i>${L[l].n}<kbd>${L[l].k}</kbd><span class="n" id="n_${l}">0</span>`;
    b.onclick = () => pick(l); bar.appendChild(b);
  });
  const cl = document.createElement('button'); cl.dataset.l = ''; cl.innerHTML = 'clear<kbd>0</kbd>';
  cl.onclick = () => pick(''); bar.appendChild(cl);
  const main = $('main');
  D.sections.forEach((sec, k) => {
    const det = document.createElement('details'); det.open = k === 0;
    const sum = document.createElement('summary'); sum.innerHTML = '<b></b>'; sum.firstChild.textContent = sec.title;
    det.appendChild(sum);
    sec.groups.forEach(g => {
      const h = document.createElement('h3'); h.textContent = (g.title ? g.title + ' ' : '') + `(${g.ids.length}) `;
      const rb = document.createElement('button'); rb.textContent = 'label the rest';
      rb.title = 'give the selected label to every unlabelled class of this section';
      rb.onclick = e => {
        e.stopPropagation();
        if (!brush) return;
        g.ids.forEach(id => { if (!st[id]) { st[id] = {l: brush, t: new Date().toISOString()}; paint(id); } });
        dirty = true; save(); counts();
      };
      h.appendChild(rb);
      const div = document.createElement('div'); div.className = 'g';
      g.ids.forEach(id => div.appendChild(tile(items[id])));
      det.append(h, div);
    });
    main.appendChild(det);
  });
  D.items.forEach(it => paint(it.id));
  pick('secondary_structure'); counts(); note();
}

function exportTSV() {
  const cols = ['dataset', 'run', 'stk', 'source', 'cls', 'pop', 'res', 'label', 'cnn_group', 'blind', 'curator', 'time'];
  const cur = $('cur').value.trim().replace(/\s+/g, ' ');
  const lines = [cols.join('\t')];
  D.items.forEach(it => {
    const s = st[it.id];
    lines.push([it.dataset, it.run, it.stk, it.source, it.cls, it.pop, it.res, s ? s.l : '', it.cnn || '',
                D.blind ? 'yes' : 'no', cur, s ? s.t : ''].join('\t'));
  });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], {type: 'text/tab-separated-values'}));
  a.download = D.fname; document.body.appendChild(a); a.click(); a.remove();
  dirty = false; msg(`exported ${D.fname}`);
}
function importTSV(file) {
  const r = new FileReader();
  r.onload = () => {
    const rows = r.result.split(/\r?\n/).filter(x => x.length);
    const h = rows.length ? rows[0].split('\t') : [], ix = k => h.indexOf(k);
    if (ix('stk') < 0 || ix('cls') < 0 || ix('label') < 0) { msg('not a label export: needs stk, cls and label columns'); return; }
    let n = 0, off = 0;
    rows.slice(1).forEach(line => {
      const f = line.split('\t'), id = f[ix('stk')] + '#' + f[ix('cls')], l = f[ix('label')];
      if (!l) return;
      if (!(id in items) || !(l in L)) { off++; return; }
      st[id] = {l: l, t: (ix('time') >= 0 && f[ix('time')]) || new Date().toISOString()}; paint(id); n++;
    });
    if (!$('cur').value && ix('curator') >= 0 && rows.length > 1) { $('cur').value = rows[1].split('\t')[ix('curator')]; }
    save(); counts(); msg(`imported ${n} labels` + (off ? `; ${off} not on this page` : ''));
  };
  r.readAsText(file);
}

document.addEventListener('keydown', e => {
  if (e.target.tagName === 'INPUT' || e.ctrlKey || e.metaKey || e.altKey) return;
  if (e.key === '0') { pick(''); return; }
  Object.keys(L).forEach(l => { if (e.key === L[l].k) pick(l); });
});
window.addEventListener('beforeunload', e => { if (dirty) { e.preventDefault(); e.returnValue = ''; } });
['ts', 'br', 'ct'].forEach(k => {
  const el = $(k);
  el.oninput = () => { document.documentElement.style.setProperty('--' + k, k === 'ts' ? el.value + 'px' : el.value); };
});
$('cur').value = (store && localStorage.getItem(CKEY)) || '';
$('cur').onchange = () => { if (store) localStorage.setItem(CKEY, $('cur').value); };
$('exp').onclick = exportTSV;
$('imp').onchange = e => { if (e.target.files.length) importTSV(e.target.files[0]); e.target.value = ''; };
build();
'''

BODY = '''<div class="top"><div class="bar"><span id="brushes" style="display:flex;flex-wrap:wrap;gap:6px"></span>
<span id="prog" class="n"></span></div>
<div class="bar"><label>Curator <input id="cur" size="12"></label>
<button id="exp">Export TSV</button><label><button onclick="document.getElementById('imp').click()">Import TSV</button>
<input id="imp" type="file" accept=".tsv,.txt" style="display:none"></label>
<label>Tile <input type="range" id="ts" min="64" max="320" step="8" value="128"></label>
<label>Brightness <input type="range" id="br" min="0.2" max="3" step="0.05" value="1"></label>
<label>Contrast <input type="range" id="ct" min="0.2" max="4" step="0.05" value="1"></label><span id="msg"></span></div></div>
'''


def page(a, items, sections):
    keys = '\n'.join(sorted(it['id'] for it in items)) + ('\nblind' if a.blind else '')
    datasets = sorted({it['dataset'] for it in items})
    fname = 'cavg_labels_{}{}.tsv'.format(re.sub(r'[^\w.-]+', '_', datasets[0]) if len(datasets) == 1 else 'mixed',
                                          '_blind' if a.blind else '')
    data = {'page': hashlib.sha1(keys.encode()).hexdigest()[:16], 'blind': a.blind, 'fname': fname, 'items': items,
            'sections': sections}
    title = 'Label class averages' + (' (blind)' if a.blind else '')
    intro = ('<p>Pick a label (button or key), then click class averages to give it; click again to clear. '
             '<b>Secondary structure</b>: helices or sheets visible. <b>Likely protein</b>: protein, no secondary '
             'structure visible. <b>Protein, ss unsure</b>: protein, but you cannot tell whether secondary structure '
             'shows. <b>Junk</b>: ice, carbon, aggregates, contamination, noise. <b>Can\'t tell</b>: protein or junk '
             'unclear. Label every class of a run, not only the ones the CNN got wrong. Labels stay in this browser; '
             'export them as TSV.' + (' Where the CNN has run, the tiles are grouped by its call and sorted by its '
                                      'probability; hover a tile for all three.' if any('cnn' in it for it in items)
                                      else '') + '</p>')
    style = STYLE % {'s': MS.SURFACE, 't1': MS.TEXT1, 't2': MS.TEXT2, 'h': MS.HAIR}
    return (f'<!doctype html><html><head><meta charset="utf-8"><title>{html.escape(title)}</title><style>{style}</style>'
            f'</head><body><h1>{html.escape(title)}</h1>{intro}{BODY}<div id="main"></div>'
            f'<script type="application/json" id="data">{json.dumps(data).replace("</", "<" + chr(92) + "/")}</script>'
            f'<script>{SCRIPT}</script></body></html>')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paths', nargs='+', help=f'projects (.simple), {TABLE} files, or folders holding tables')
    ap.add_argument('--out', required=True, help='the labelling page (HTML)')
    ap.add_argument('--dataset', help='the data\'s name in the labels (default: the project\'s name)')
    ap.add_argument('--blind', action='store_true', help='leave the CNN\'s calls out of the page')
    ap.add_argument('--simple', help='SIMPLE build or install folder holding bin/simple_exec (default $SIMPLE_PATH)')
    a = ap.parse_args()
    sp = MS.find_simple(a.simple)
    env = dict(os.environ, SIMPLE_PATH=sp)
    env.setdefault('SIMPLE_QSYS', 'local')
    items, sections = collect(a, os.path.join(sp, 'bin', 'simple_exec'), env)
    if not items:
        raise SystemExit('no classes with particles found')
    with open(a.out, 'w') as f:
        f.write(page(a, items, sections))
    print(f'{len(items)} classes in {len(sections)} run{"" if len(sections) == 1 else "s"}: {os.path.abspath(a.out)}')


if __name__ == '__main__':
    sys.exit(main())
