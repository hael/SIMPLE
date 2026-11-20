#!/usr/bin/env python3
"""Search EMDB for map candidates: single-particle entries with both half maps, a fitted model, a resolution better
than --maxres, a map no larger than --maxbox and, with --min-kda/--max-kda, an overall molecular weight in that range;
printed one per line (id, resolution, weight, box, voxel, field, model, title).

Usage: python -I emdb_search.py --query "MscS AND nanodisc" [--maxres 3.6] [--maxbox 480] [--rows 60] [--pages 1]
       [--min-kda 30 --max-kda 150]
EMDB returns at most 100 entries per page; --pages N reads up to N pages.
--query may be any EMDB search (Solr) expression, e.g. 'structure_determination_method:"singleParticle"'.
Standard library only; queries https://www.ebi.ac.uk/emdb/api/search.
"""
import argparse
import json
import sys
import urllib.parse
import urllib.request


def weight_kda(e):
    """The sample's overall theoretical weight, or the sum over its macromolecules above 2 kDa (no water, ions)."""
    unit = {'MDa': 1000.0, 'kDa': 1.0, 'Da': 0.001}
    s = e.get('sample') or {}
    for x in (s.get('supramolecule_list') or {}).get('supramolecule') or []:
        w = (x.get('molecular_weight') or {}).get('theoretical') or {}
        if w.get('valueOf_'):
            return float(w['valueOf_']) * unit.get(w.get('units'), 1.0)
    tot = 0.0
    for m in (s.get('macromolecule_list') or {}).get('macromolecule') or []:
        w = (m.get('molecular_weight') or {}).get('theoretical') or {}
        if w.get('valueOf_'):
            kda = float(w['valueOf_']) * unit.get(w.get('units'), 1.0)
            if kda > 2:
                tot += kda * float(m.get('number_of_copies') or 1)
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--query', required=True)
    ap.add_argument('--maxres', type=float, default=3.6)
    ap.add_argument('--maxbox', type=int, default=480)
    ap.add_argument('--rows', type=int, default=60, help='entries per page (EMDB allows 100) {60}')
    ap.add_argument('--pages', type=int, default=1, help='pages to read {1}')
    ap.add_argument('--min-kda', type=float, default=0.0, help='least overall molecular weight (kDa)')
    ap.add_argument('--max-kda', type=float, default=0.0, help='largest overall molecular weight (kDa; 0: no limit)')
    a = ap.parse_args()
    q = f'({a.query}) AND resolution:[0 TO {a.maxres}]'
    if a.min_kda or a.max_kda:
        q += f' AND overall_molecular_weight:[{a.min_kda / 1000} TO {(a.max_kda or 1e6) / 1000}]'
    data = []
    for page in range(1, a.pages + 1):
        url = 'https://www.ebi.ac.uk/emdb/api/search/' + urllib.parse.quote(q) + f'?rows={a.rows}&page={page}'
        with urllib.request.urlopen(url, timeout=300) as r:
            got = json.load(r)
        data += got
        if len(got) < a.rows:
            break
    out = []
    for e in data:
        try:
            sd = e['structure_determination_list']['structure_determination'][0]
            method = sd.get('method', '')
            res = float(sd['image_processing'][0]['final_reconstruction']['resolution']['valueOf_'])
            halves = len(((e.get('interpretation') or {}).get('half_map_list') or {}).get('half_map', []) or [])
            pdbs = [p.get('pdb_id') for p in ((e.get('crossreferences') or {}).get('pdb_list') or {}).get('pdb_reference', [])]
            m = e['map']
            box = int(m['dimensions']['col'])
            apix = float(m['pixel_spacing']['x']['valueOf_'])
        except Exception:
            continue
        if method != 'singleParticle' or halves < 2 or not pdbs or res >= a.maxres or box > a.maxbox:
            continue
        kda = weight_kda(e)
        if (a.min_kda and kda < a.min_kda) or (a.max_kda and kda > a.max_kda):
            continue
        out.append((res, e['emdb_id'], kda, box, apix, pdbs[0], e['admin']['title'][:90]))
    for res, eid, kda, box, apix, pdb, title in sorted(out):
        print(f'{eid:10s} {res:4.2f} A  {kda:5.0f} kDa  box {box:3d} x {apix:.3f} A = {box * apix:5.0f} A  model {pdb}  | {title}')
    print(f'{len(out)} of {len(data)} entries found pass', file=sys.stderr)


if __name__ == '__main__':
    main()
