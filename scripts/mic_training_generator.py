#!/usr/bin/env python3
"""solve2D on the particles of the first N micrographs of a project, for several N, the class-average CNN
(SIMPLE's cavg_cnn_predict) on each result, and one report over the series.

  python3 mic_training_generator.py --project proj.simple --mskdiam 175 --out DIR [--nmics 10,20,50,100,200] [--nthr 16]
                        [--simple SIMPLE_BUILD] [--weights W.bin] [--ncls 0] [--ptcls-per-class 200]

Per N, in DIR/mics_<N>/:
  1. selection oritype=mic with the first N micrographs that are selected and hold particles flagged 1
     (mic_flags.txt), so their stacks' particles stay selected and the others get state 0 -> 1_selection
  2. solve2D on that project -> 2_solve2D; ncls is --ncls, or the selected particles / --ptcls-per-class kept
     within [--min-ncls, --max-ncls]
  3. cavg_cnn_predict on the result -> 3_cavg_cnn_predict/cavg_cnn_predictions.txt
The input project is never written: each program works on a copy in its own execution folder. A step that finished
is not run again, so the series resumes where it stopped. The project needs one particle stack per micrograph (as
the stream writes them); otherwise selection cannot carry micrograph states to the particles and the series stops.
DIR/report.html: per N the micrographs, particles, classes, the CNN's groups by classes and particles, the share kept
as protein, the run times, then the class averages grouped by the CNN's call. DIR/series.json holds the numbers.
Standard library only (Python 3.7 or newer); SIMPLE from --simple, else $SIMPLE_PATH, else the install or the
repository's build/ this script sits in.
"""
import argparse
import array
import base64
import glob
import html
import json
import os
import re
import struct
import subprocess
import sys
import time
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
GROUPS = ('secondary_structure', 'likely_protein', 'junk', 'empty')
TITLES = {'secondary_structure': 'secondary structure', 'likely_protein': 'likely protein', 'junk': 'junk',
          'empty': 'empty'}
PROB = {'secondary_structure': 'p_ss', 'likely_protein': 'p_likely', 'junk': 'p_junk'}
SURFACE, TEXT1, TEXT2, HAIR = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'


def say(msg):
    print('[series {}] {}'.format(time.strftime('%H:%M:%S'), msg))
    sys.stdout.flush()


def find_simple(given):
    """--simple, else $SIMPLE_PATH, else the SIMPLE install holding this script, else the repository's build/."""
    for sp in (given, os.environ.get('SIMPLE_PATH'), os.path.join(HERE, '..'), os.path.join(HERE, '..', 'build')):
        if sp and os.access(os.path.join(sp, 'bin', 'simple_exec'), os.X_OK):
            return os.path.abspath(sp)
        if sp and sp in (given, os.environ.get('SIMPLE_PATH')):
            raise SystemExit('no {}'.format(os.path.join(sp, 'bin', 'simple_exec')))
    raise SystemExit('SIMPLE not found: give --simple, the SIMPLE build or install folder holding bin/simple_exec')


# ------------------------------------------------------------------ the runs
def records(exe, env, projfile, oritype, cwd):
    """The key=value records of a project segment, as dicts of strings."""
    r = subprocess.run([exe, 'prg=print_project_field', 'projfile=' + projfile, 'oritype=' + oritype], cwd=cwd, env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    return [dict(re.findall(r'(\w+)=(\S+)', l)) for l in r.stdout.splitlines() if '=' in l]


def run(exe, env, args, cwd, log, stop):
    """A SIMPLE program in cwd, its output in log; True when the log holds its normal stop (SIMPLE can exit 0 on
    errors, so the exit code alone is not trusted)."""
    t0 = time.time()
    with open(os.path.join(cwd, log), 'w') as f:
        rc = subprocess.run([exe] + args, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT).returncode
    with open(os.path.join(cwd, log), errors='replace') as f:
        ok = rc == 0 and stop in f.read()
    return ok, time.time() - t0


def newest(cwd, pattern):
    """The highest-numbered execution folder matching pattern."""
    hits = [p for p in glob.glob(os.path.join(cwd, pattern)) if os.path.basename(p).split('_')[0].isdigit()]
    hits.sort(key=lambda p: int(os.path.basename(p).split('_')[0]))
    return hits[-1] if hits else None


def save(rec, path):
    with open(path, 'w') as f:
        json.dump(rec, f, indent=1)


def one_n(a, exe, env, n, mics_ok, proj_name, out):
    """The three steps for the first n eligible micrographs; returns the run's record."""
    d = os.path.join(out, 'mics_{}'.format(n))
    os.makedirs(d, exist_ok=True)
    rec_p = os.path.join(d, 'record.json')
    rec = {'nmics': n}
    if os.path.exists(rec_p):
        with open(rec_p) as f:
            rec = json.load(f)
    if rec.get('cnn_ok'):
        say('{} micrographs: done before'.format(n))
        return rec
    # 1. the selection
    if not rec.get('selection_ok'):
        flags, k = [], 0
        for ok in mics_ok:
            flags.append(1 if ok and k < n else 0)
            k += ok
        with open(os.path.join(d, 'mic_flags.txt'), 'w') as f:
            f.write('\n'.join(str(x) for x in flags) + '\n')
        say('{} micrographs: selection'.format(n))
        ok, _ = run(exe, env, ['prg=selection', 'projfile=' + os.path.abspath(a.project), 'oritype=mic',
                               'infile=mic_flags.txt'], d, 'selection.log', 'SELECTION NORMAL STOP')
        sel = newest(d, '*_selection')
        if not ok or not sel:
            rec['error'] = 'selection failed; see selection.log'
            save(rec, rec_p)
            return rec
        rec['selection_proj'] = os.path.join(sel, proj_name)
        ptcls = records(exe, env, rec['selection_proj'], 'ptcl2D', d)
        rec['nptcls'] = sum(1 for p in ptcls if float(p.get('state', 1)) > 0)
        rec['selection_ok'] = True
        save(rec, rec_p)
    # 2. solve2D
    if not rec.get('solve2D_ok'):
        ncls = a.ncls if a.ncls else int(min(a.max_ncls, max(a.min_ncls, round(rec['nptcls'] / a.ptcls_per_class))))
        rec['ncls'] = ncls
        say('{} micrographs: solve2D on {} particles, {} classes'.format(n, rec['nptcls'], ncls))
        ok, t = run(exe, env, ['prg=solve2D', 'projfile=' + rec['selection_proj'], 'ncls={}'.format(ncls),
                               'mskdiam={}'.format(a.mskdiam), 'nthr={}'.format(a.nthr)], d, 'solve2D.log',
                    'SOLVE2D NORMAL STOP')
        s2 = newest(d, '*_solve2D')
        rec['solve2D_s'] = round(t, 1)
        if not ok or not s2:
            rec['error'] = 'solve2D failed; see solve2D.log'
            save(rec, rec_p)
            return rec
        rec['solve2D_proj'] = os.path.join(s2, proj_name)
        rec['solve2D_ok'] = True
        save(rec, rec_p)
    # 3. the class-average CNN
    say('{} micrographs: cavg_cnn_predict'.format(n))
    ok, t = run(exe, env, ['prg=cavg_cnn_predict', 'projfile=' + rec['solve2D_proj'],
                           'infile=' + os.path.abspath(a.weights), 'nthr={}'.format(a.nthr)], d,
                'cavg_cnn_predict.log', 'CAVG_CNN_PREDICT NORMAL STOP')
    cn = newest(d, '*_cavg_cnn_predict')
    rec['cnn_s'] = round(t, 1)
    if not ok or not cn:
        rec['error'] = 'cavg_cnn_predict failed; see cavg_cnn_predict.log'
        save(rec, rec_p)
        return rec
    rec['table'] = os.path.join(cn, 'cavg_cnn_predictions.txt')
    stk = [o.get('stk') for o in records(exe, env, rec['solve2D_proj'], 'out', d) if o.get('imgkind') == 'cavg']
    rec['cavgs'] = os.path.join(os.path.dirname(rec['solve2D_proj']), stk[0]) if stk else None
    rec['cnn_ok'] = True
    rec.pop('error', None)
    save(rec, rec_p)
    return rec


# ------------------------------------------------------------------ the report
def read_mrc(path):
    """The images of a 32-bit MRC stack, each a flat array of n*n floats, and n."""
    with open(path, 'rb') as f:
        h = f.read(1024)
        nx, ny, nz, mode = struct.unpack('<4i', h[:16])
        if mode != 2 or nx != ny:
            raise ValueError('{}: not a stack of square 32-bit images'.format(path))
        f.seek(1024 + struct.unpack('<i', h[92:96])[0])
        a = array.array('f')
        a.fromfile(f, nx * ny * nz)
    if sys.byteorder == 'big':
        a.byteswap()
    return [a[k * nx * ny:(k + 1) * nx * ny] for k in range(nz)], nx


def png(img, n, size=96):
    """A size x size greyscale PNG of an n x n image, nearest samples, its minimum black and its maximum white."""
    idx = [int(i * n / size) for i in range(size)]
    vals = [img[j * n + i] for j in idx for i in idx]
    lo = min(vals)
    span = (max(vals) - lo) or 1.0
    g = bytes(min(255, max(0, int((v - lo) / span * 255))) for v in vals)
    raw = b''.join(b'\x00' + g[j * size:(j + 1) * size] for j in range(size))

    def chunk(tag, data):
        c = tag + data
        return struct.pack('>I', len(data)) + c + struct.pack('>I', zlib.crc32(c) & 0xffffffff)
    return (b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', size, size, 8, 0, 0, 0, 0))
            + chunk(b'IDAT', zlib.compress(raw, 6)) + chunk(b'IEND', b''))


def tile(img, n, cap, tip):
    src = 'data:image/png;base64,' + base64.b64encode(png(img, n)).decode()
    return ('<div class="t"><img class="cv" src="{}" width="96" height="96" title="{}" style="display:block;'
            'border-radius:3px"><div class="cap">{}</div></div>').format(src, html.escape(tip), cap)


def table(header, rows):
    return ('<table><tr>' + ''.join('<th>{}</th>'.format(html.escape(h)) for h in header) + '</tr>'
            + ''.join('<tr>' + ''.join('<td>{}</td>'.format(html.escape(str(v))) for v in r) + '</tr>' for r in rows)
            + '</table>')


VIEW_BAR = ('<div class="bar"><label>Brightness <input type="range" id="br" min="0.2" max="3" step="0.05" value="1">'
            '<span id="brv">1.00</span></label><label>Contrast <input type="range" id="ct" min="0.2" max="4" step="0.05" '
            'value="1"><span id="ctv">1.00</span></label><button id="rs">Reset</button></div>'
            '<script>(function(){var r=document.documentElement.style;'
            'function u(k,v){r.setProperty("--"+k,v);document.getElementById(k+"v").textContent=(+v).toFixed(2);}'
            '["br","ct"].forEach(function(k){var e=document.getElementById(k);e.oninput=function(){u(k,e.value);};});'
            'document.getElementById("rs").onclick=function(){["br","ct"].forEach(function(k){'
            'document.getElementById(k).value=1;u(k,1);});};})();</script>')


def page(title, body):
    css = ('body{{background:{s};color:{t1};font:14px/1.45 system-ui,sans-serif;margin:24px;max-width:1200px}}'
           'h1{{font-size:20px}}h2{{font-size:16px;margin-top:28px}}h3{{font-size:14px}}p,.cap{{color:{t2}}}'
           'table{{border-collapse:collapse;margin:8px 0}}td,th{{border-bottom:1px solid {h};padding:3px 10px;'
           'text-align:right;font-variant-numeric:tabular-nums}}th{{color:{t2};font-weight:600}}'
           'td:first-child,th:first-child{{text-align:left}}.g{{display:flex;flex-wrap:wrap;gap:10px}}'
           '.t{{width:96px}}.cap{{font-size:11px;line-height:1.25;margin-top:2px}}summary{{cursor:pointer;margin:10px 0}}'
           'img.cv{{filter:brightness(var(--br,1)) contrast(var(--ct,1))}}'
           '.bar{{position:sticky;top:0;z-index:2;background:{s};border-bottom:1px solid {h};padding:8px 0;'
           'display:flex;gap:24px;align-items:center;color:{t2}}}.bar label{{display:flex;gap:8px;align-items:center}}'
           '.bar span{{min-width:2.6em;font-variant-numeric:tabular-nums;color:{t1}}}').format(
        s=SURFACE, t1=TEXT1, t2=TEXT2, h=HAIR)
    return ('<!doctype html><html><head><meta charset="utf-8"><title>{0}</title><style>{1}</style></head><body>'
            '<h1>{0}</h1>{2}{3}</body></html>').format(html.escape(title), css, VIEW_BAR, body)


def summarise(rec):
    """The run's table rows and, per group, its classes and particles."""
    with open(rec['table']) as f:
        rows = [dict(re.findall(r'(\w+)=(\S+)', l)) for l in f if '=' in l]
    g = {k: {'classes': 0, 'particles': 0} for k in GROUPS}
    for r in rows:
        g[r['group']]['classes'] += 1
        g[r['group']]['particles'] += int(r['pop'])
    return rows, g


def report(a, recs, out):
    head = ['micrographs', 'particles', 'classes (with particles)'] + \
        ['{}: classes / particles'.format(TITLES[k]) for k in GROUPS[:3]] + ['kept as protein', 'solve2D', 'CNN']
    rows, sections = [], ''
    done = [r for r in recs if r.get('cnn_ok')]
    for rec in recs:
        n = rec['nmics']
        if not rec.get('cnn_ok'):
            rows.append([n, rec.get('nptcls', ''), '', '', '', '', rec.get('error', 'not run'), '', ''])
            continue
        trows, g = summarise(rec)
        tot = sum(g[k]['particles'] for k in GROUPS)
        kept = g['secondary_structure']['particles'] + g['likely_protein']['particles']
        rows.append([n, rec['nptcls'], '{} ({})'.format(rec['ncls'], rec['ncls'] - g['empty']['classes'])] +
                    ['{} / {}'.format(g[k]['classes'], g[k]['particles']) for k in GROUPS[:3]] +
                    ['{} ({:.0f}%)'.format(kept, 100 * kept / max(tot, 1)), '{:.1f} min'.format(rec['solve2D_s'] / 60),
                     '{:.1f} s'.format(rec['cnn_s'])])
        if not rec.get('cavgs') or not os.path.exists(rec['cavgs']):
            continue
        imgs, box = read_mrc(rec['cavgs'])
        body = ''
        for k in GROUPS[:3]:
            sel = sorted((r for r in trows if r['group'] == k), key=lambda r: -float(r[PROB[k]]))
            tiles = [tile(imgs[int(r['class']) - 1], box,
                          'c{} · p {:.2f} · {} ptcls<br>{:.1f} Å'.format(r['class'], float(r[PROB[k]]), r['pop'],
                                                                         float(r['res'])),
                          'class {}: junk {:.3f}, ss {:.3f}, likely {:.3f}'.format(
                              r['class'], float(r['p_junk']), float(r['p_ss']), float(r['p_likely']))) for r in sel]
            body += '<h3>{}: {} class{}, {} particles</h3><div class="g">{}</div>'.format(
                TITLES[k], len(sel), '' if len(sel) == 1 else 'es', g[k]['particles'], ''.join(tiles))
        sections += '<details{}><summary><b>{} micrographs</b>: {} particles, {} classes</summary>{}</details>'.format(
            ' open' if rec is done[-1] else '', n, rec['nptcls'], rec['ncls'], body)
    classes = str(a.ncls) if a.ncls else 'particles / {:g} within [{}, {}]'.format(a.ptcls_per_class, a.min_ncls, a.max_ncls)
    intro = ('<p>{}; mask {:g} Å; {} threads; classes {}; weights {}; {}.</p><p>Each N is the first N micrographs of '
             'the project that are selected and hold particles. Kept as protein: the particles of the classes the CNN '
             'calls secondary structure or likely protein. Class captions: class, the probability of its group, '
             'particles, resolution; hover for all three probabilities.</p>').format(
        html.escape(os.path.abspath(a.project)), a.mskdiam, a.nthr, html.escape(classes),
        html.escape(os.path.abspath(a.weights)), time.strftime('%Y-%m-%d %H:%M'))
    with open(os.path.join(out, 'report.html'), 'w') as f:
        f.write(page('Class-average CNN over a micrograph series',
                     intro + table(head, rows) + '<h2>Class averages</h2>' + sections))
    say('report: ' + os.path.join(out, 'report.html'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--project', required=True, help='a SIMPLE project with micrographs and their particle stacks')
    ap.add_argument('--mskdiam', type=float, required=True, help='mask diameter for solve2D (A)')
    ap.add_argument('--out', required=True, help='output folder')
    ap.add_argument('--nmics', default='10,20,50,100,200', help='micrograph counts (default %(default)s)')
    ap.add_argument('--nthr', type=int, default=16, help='threads per SIMPLE program (default %(default)s)')
    ap.add_argument('--simple', help='SIMPLE build or install folder holding bin/simple_exec (default $SIMPLE_PATH)')
    ap.add_argument('--weights', default=os.path.join(HERE, 'sscnn', 'weights', 'cavg_cnn_v3.bin'),
                    help='the CNN weights for cavg_cnn_predict (default sscnn/weights/cavg_cnn_v3.bin beside this script)')
    ap.add_argument('--ncls', type=int, default=0, help='solve2D classes; 0 (default): from the particle count')
    ap.add_argument('--ptcls-per-class', type=float, default=200., help='(default %(default)s)')
    ap.add_argument('--min-ncls', type=int, default=10, help='(default %(default)s)')
    ap.add_argument('--max-ncls', type=int, default=100, help='(default %(default)s)')
    a = ap.parse_args()
    sp = find_simple(a.simple)
    exe = os.path.join(sp, 'bin', 'simple_exec')
    env = dict(os.environ, SIMPLE_PATH=sp)
    env.setdefault('SIMPLE_QSYS', 'local')
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    if not os.path.exists(a.project):
        raise SystemExit('no project {}'.format(a.project))
    if not os.path.exists(a.weights):
        raise SystemExit('no weights file {}: give --weights'.format(a.weights))
    say('SIMPLE {}'.format(sp))
    mics = records(exe, env, os.path.abspath(a.project), 'mic', out)
    stks = records(exe, env, os.path.abspath(a.project), 'stk', out)
    if not mics:
        raise SystemExit('the project has no micrographs')
    if len(stks) != len(mics):
        raise SystemExit('{} micrographs but {} particle stacks: selection needs one stack per micrograph'.format(
            len(mics), len(stks)))
    # selection carries micrograph i's state to stack i's particles, so stack i's count is micrograph i's
    mics_ok = [float(m.get('state', 1)) > 0 and float(s.get('state', 1)) > 0 and float(s.get('nptcls', 0)) > 0
               for m, s in zip(mics, stks)]
    say('{} micrographs, {} selected with particles'.format(len(mics), sum(mics_ok)))
    proj_name = os.path.basename(a.project)
    recs = []
    for n in sorted(int(x) for x in a.nmics.split(',')):
        if n > sum(mics_ok):
            say('{} micrographs: skipped, the project has {}'.format(n, sum(mics_ok)))
            recs.append({'nmics': n, 'error': 'the project has only {} micrographs with particles'.format(sum(mics_ok))})
        else:
            recs.append(one_n(a, exe, env, n, mics_ok, proj_name, out))
        save(recs, os.path.join(out, 'series.json'))
        report(a, recs, out)                                  # after every N, so a long series can be read as it goes
    say('done')


if __name__ == '__main__':
    main()
