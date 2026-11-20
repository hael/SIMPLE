#!/usr/bin/env python3
"""An HTML report of simulated runs, images embedded: the input maps (projections and central slices), a table of
per-run numbers, and per run its most populated protein classes, each with its helix truth as it survives in the class
average, the band-passed helix and strand truth (the training targets), then its junk classes.

Band-passed truth is shown on one grey scale for the whole report (the median spread of the first run's helix targets),
so that weak targets look weak.

Usage: python -I report.py --title T --out report.html --run LABEL=RUN_DIR[=MAP_DIR] [--run ...] [--nclasses 8] [--njunk 6]
       [--intro TEXT]
"""
import argparse
import base64
import csv
import html
import json
import os
import struct
import sys
import zlib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ss_trainer as T   # noqa: E402
from sim_truth import load_truth   # noqa: E402

SURFACE, TEXT1, TEXT2, HAIR = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'


def png_bytes(g):
    h, w = g.shape
    raw = b''.join(b'\x00' + g[j].tobytes() for j in range(h))

    def ch(t, d):
        return struct.pack('>I', len(d)) + t + d + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
    return b'\x89PNG\r\n\x1a\n' + ch(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 0, 0, 0, 0)) + ch(b'IDAT', zlib.compress(raw, 6)) + ch(b'IEND', b'')


def to_uri(g):
    return 'data:image/png;base64,' + base64.b64encode(png_bytes(g)).decode()


def auto(img, size):
    """Grey over +-3 sd of the image itself."""
    return T.png(np.asarray(img, np.float64), size)


def fixed(img, centre, sd, size):
    """Grey over centre +-3 sd, a scale shared by several images."""
    n = img.shape[0]
    idx = (np.arange(size) * n / size).astype(int)
    im = img[np.ix_(idx, idx)]
    return np.clip((im - centre + 3 * sd) / (6 * sd) * 255, 0, 255).astype(np.uint8)


def uri_auto(img, size):
    return 'data:image/png;base64,' + base64.b64encode(auto(img, size)).decode()


def read_vol(path):
    with open(path, 'rb') as f:
        h = f.read(1024)
        n = struct.unpack('<i', h[:4])[0]
        return np.fromfile(f, '<f4', n ** 3).reshape(n, n, n).astype(np.float64)


def img_tag(uri, w, title=''):
    return f'<img src="{uri}" width="{w}" height="{w}" title="{html.escape(title)}" style="display:block;border-radius:3px">'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--title', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--run', action='append', required=True)
    ap.add_argument('--nclasses', type=int, default=8)
    ap.add_argument('--njunk', type=int, default=6)
    ap.add_argument('--intro', default='')
    a = ap.parse_args()
    runs = []
    for spec in a.run:
        parts = spec.split('=')
        runs.append((parts[0], parts[1], parts[2] if len(parts) > 2 else ''))
    # the shared scale of the band-passed targets
    first = runs[0][1]
    rows0 = [r for r in csv.DictReader(open(os.path.join(first, 'classes.tsv')), delimiter='\t')
             if int(r['pop']) > 0 and int(r['njunk']) <= 0.5 * int(r['pop'])]
    b0 = load_truth(first, 'helix_band')
    sd_band = float(np.median([b0[int(r['class']) - 1].std() for r in rows0]))
    # maps panel
    maps_html = ''
    mp = [(lab, md) for lab, _, md in runs if md and os.path.exists(os.path.join(md, 'full.mrc'))]
    seen = set()
    cells = []
    for lab, md in mp:
        if md in seen:
            continue
        seen.add(md)
        v = read_vol(os.path.join(md, 'full.mrc'))
        n = v.shape[0]
        cells.append(f'<figure style="margin:0"><div style="display:flex;gap:4px">{img_tag(uri_auto(v.sum(0), 150), 150, "projection along z")}'
                     f'{img_tag(uri_auto(v.sum(2), 150), 150, "projection along x")}{img_tag(uri_auto(v[:, :, n // 2], 150), 150, "central slice")}</div>'
                     f'<figcaption style="font-size:11px;color:{TEXT2};margin-top:3px">{html.escape(lab)}: projection down z, down x, central slice</figcaption></figure>')
    if cells:
        maps_html = (f'<h2>Input maps</h2><div style="display:flex;flex-wrap:wrap;gap:18px">{"".join(cells)}</div>')
    # numbers
    trs = []
    sections = []
    for lab, rd, _ in runs:
        info = json.load(open(os.path.join(rd, 'run.json')))
        rows = list(csv.DictReader(open(os.path.join(rd, 'classes.tsv')), delimiter='\t'))
        prot = [r for r in rows if int(r['pop']) > 0 and int(r['njunk']) <= 0.5 * int(r['pop'])]
        junk = [r for r in rows if int(r['pop']) > 0 and int(r['njunk']) > 0.5 * int(r['pop'])]
        cav, smpd, _, _ = T.read_mrc(os.path.join(rd, 'cavgs.mrc'))
        ev = T.read_mrc(os.path.join(rd, 'cavgs_even.mrc'))[0]
        od = T.read_mrc(os.path.join(rd, 'cavgs_odd.mrc'))[0]
        zero = np.zeros_like(cav)
        tf, th, hb, sb = [x if x is not None else zero for x in
                          (load_truth(rd, k) for k in ('full_on_cavg', 'helix_on_cavg', 'helix_band', 'strand_band'))]
        M = T.Measurer(cav.shape[1], smpd, info['mskdiam'])
        feats = np.array([M.measure(cav[int(r['class']) - 1], ev[int(r['class']) - 1], od[int(r['class']) - 1]) for r in prot])
        med = lambda k: float(np.median([float(r[k]) for r in prot])) if prot else float('nan')
        band_share = float(np.median([hb[int(r['class']) - 1].std() / max(cav[int(r['class']) - 1].std(), 1e-30) for r in prot]))
        jk = json.load(open(os.path.join(rd, 'junk.json')))
        trs.append(f'<tr><td>{html.escape(lab)}</td><td>{info["snr"]:g}</td><td>{len(prot)}</td><td>{len(junk)} '
                   f'({sum(int(r["njunk"]) for r in junk)} of {jk["njunk"]} junk)</td><td>{med("res_A"):.1f}</td>'
                   f'<td>{med("cnr_helix_12_6"):.2f}</td><td>{sum(float(r["cnr_helix_12_6"]) >= 1 for r in prot)}</td>'
                   f'<td>{med("cnr_strand_12_6"):.2f}</td><td>{band_share:.3f}</td><td>{np.median(feats[:, 1]):.1f}</td>'
                   f'<td>{np.median(feats[:, 0]):.2f}</td></tr>')
        prot.sort(key=lambda r: -int(r['pop']))
        junk.sort(key=lambda r: -int(r['pop']))
        pc = []
        for r in prot[:a.nclasses]:
            c = int(r['class']) - 1
            s_full, m_full = tf[c].std() or 1.0, tf[c].mean()
            pc.append(
                f'<figure style="margin:0;width:330px"><div style="display:flex;gap:3px;align-items:flex-start">'
                f'{img_tag(uri_auto(cav[c], 120), 120, "class average")}'
                f'{img_tag(to_uri(fixed(th[c], m_full, s_full, 66)), 66, "helix truth in the class average")}'
                f'{img_tag(to_uri(fixed(hb[c], 0.0, sd_band, 66)), 66, "band-passed helix truth (target)")}'
                f'{img_tag(to_uri(fixed(sb[c], 0.0, sd_band, 66)), 66, "band-passed strand truth (target)")}</div>'
                f'<figcaption style="font-size:10px;color:{TEXT2};margin-top:2px">class {c + 1} · n={r["pop"]} · {float(r["res_A"]):.1f} Å · '
                f'helix cnr {float(r["cnr_helix_12_6"]):.2f} · strand cnr {float(r["cnr_strand_12_6"]):.2f}</figcaption></figure>')
        jc = ''.join(
            f'<figure style="margin:0;width:110px">{img_tag(uri_auto(cav[int(r["class"]) - 1], 104), 104)}'
            f'<figcaption style="font-size:10px;color:{TEXT2};margin-top:2px">class {r["class"]} · n={r["pop"]} · '
            f'{100 * int(r["njunk"]) / int(r["pop"]):.0f}% junk</figcaption></figure>' for r in junk[:a.njunk])
        sections.append(f'<h2>{html.escape(lab)} <span style="font-weight:400;color:{TEXT2}">· particle SNR {info["snr"]:g} · '
                        f'{info["nptcls"]} protein + {info["njunk"]} junk particles · {info["seconds"]} s</span></h2>'
                        f'<h3>Most populated protein classes</h3><div style="display:flex;flex-wrap:wrap;gap:12px">{"".join(pc)}</div>'
                        f'<h3>Most populated junk classes</h3><div style="display:flex;flex-wrap:wrap;gap:10px">{jc}</div>')
    table = (f'<table><tr><th>run</th><th>particle SNR</th><th>protein classes</th><th>junk classes</th>'
             f'<th>median resolution (Å)</th><th>helix cnr 12–6 Å (median)</th><th>classes with helix cnr ≥ 1</th>'
             f'<th>strand cnr 12–6 Å (median)</th><th>band-passed helix target / class average (sd)</th>'
             f'<th>ridge measure (trainer)</th><th>half FRC 12–6 Å (trainer)</th></tr>{"".join(trs)}</table>')
    style = (f'body{{background:{SURFACE};color:{TEXT1};font:13px system-ui,sans-serif;margin:24px;max-width:1500px}}'
             f'h1{{font-size:19px;margin:0 0 6px}}h2{{font-size:15px;margin:26px 0 6px;border-top:1px solid {HAIR};padding-top:14px}}'
             f'h3{{font-size:12px;margin:12px 0 6px;color:{TEXT2};font-weight:600}}'
             f'table{{border-collapse:collapse;font-size:12px}}td,th{{border-bottom:1px solid {HAIR};padding:5px 9px;text-align:right}}'
             f'td:first-child,th:first-child{{text-align:left}}th{{font-weight:600;color:{TEXT2};vertical-align:bottom}}'
             f'p{{max-width:1000px;line-height:1.45}}')
    legend = (f'<p style="color:{TEXT2};font-size:12px">Each protein class: the class average; its helix truth as it survives in '
              f'the class average (on the scale of the class\'s noise-free average); its band-passed helix truth and band-passed '
              f'strand truth, the training targets (finer than 12–15 Å; one grey scale for the whole report, mid-grey = zero). '
              f'cnr: the truth\'s contrast against the class\'s noise at 12–6 Å. Hover an image for its name.</p>')
    page = (f'<!doctype html><meta charset="utf-8"><title>{html.escape(a.title)}</title><style>{style}</style>'
            f'<h1>{html.escape(a.title)}</h1>{a.intro}{maps_html}<h2>Numbers</h2>{table}{legend}{"".join(sections)}')
    open(a.out, 'w').write(page)
    print('wrote', a.out, f'({os.path.getsize(a.out) / 1e6:.1f} MB)')


if __name__ == '__main__':
    main()
