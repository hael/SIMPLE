#!/usr/bin/env python3
"""test_review_dossier.py -- review dossiers and the test inventory.

Phase 0 and section 9 of doc/refactoring_notes/uniform_test_environment_refactoring.md.

Reads the tree (production/tests -- gone since the utils review, the standalone
programs are retired --, src/main/commanders/test, the test UI and
routers, CI, scripts, doc) and optional timing runs from
scripts/test_timing_run.sh, and writes

  * one dossier per test identity (Markdown) under --dossiers, and
  * the inventory skeleton, one table per area, at --inventory.

An identity is a canonical test name; it may have two routes (a standalone
program and a simple_test_exec case). Everything here is static inspection;
timings come only from a --timing file. Nothing is compiled or run.

usage:
  scripts/test_review_dossier.py [--timing build_test_runs/<label>/results.tsv ...]
                                 [--inventory doc/refactoring_notes/test_inventory.md]
                                 [--dossiers build_test_runs/dossiers]
                                 [--coverage-after INVENTORY]   # report coverage lost by verdicts

--coverage-after reads the verdict column of an inventory that reviewers have
filled in, drops every identity marked delete / merge into / retire, and lists
the production procedures and modules that no remaining test exercises
(section 9.5 of the plan).
"""
import argparse
import collections
import glob
import itertools
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RE_USE = re.compile(r'^\s*use\s+(\w+)', re.I | re.M)
RE_TBP = re.compile(r'\b(\w+)(?:\([^()]*\))?\s*%\s*(\w+)\s*\(', re.I)   # obj%meth( and arr(i)%meth(
RE_TBP_BARE = re.compile(r'^\s*call\s+(\w+)(?:\([^()]*\))?\s*%\s*(\w+)\s*$', re.I | re.M)   # call obj%meth with no argument list
RE_CALL = re.compile(r'^\s*call\s+(\w+)\s*(?:\((?![^()]*\)\s*%)|$)', re.I | re.M)   # call proc( or bare call proc; not call arr(i)%meth(
RE_SUB = re.compile(r'^\s*(?:recursive\s+)?subroutine\s+(\w+)\s*\(.*?^\s*end\s+subroutine\s+\1', re.I | re.M | re.S)
NOT_PRODUCTION = {'simple_core_module_api', 'simple_commanders_api', 'simple_test_exec_api', 'simple_test_utils',
                  'simple_cmdline', 'simple_parameters', 'simple_defs', 'simple_string', 'simple_syslib',
                  'simple_fileio', 'simple_timer', 'simple_jiffys', 'simple_error', 'simple_defs_environment',
                  'simple_local_flags', 'iso_c_binding', 'iso_fortran_env', 'omp_lib', 'simple_string_utils'}
NOT_OBJECT = {'self', 'cline', 'params', 'p', 'cline_', 'spproj', 'build', 'b'}   # objects whose calls say little about coverage
PLATFORM = [('coarray', r'\bcoarray|\bco_(sum|max|min|broadcast)\b|\bsync\s+(all|images)\b|\bthis_image\s*\(|\bnum_images\s*\(|codimension|\w\[\s*\w+\s*\]\s*[=%]'), ('mpi', r'\bmpi_'),
            ('offload', r'omp target|openmp_offload|USE_OPENMP_OFFLOAD'), ('openacc', r'!\$acc|openacc'),
            ('cuda', r'\bcuda|flex_gpu'), ('socket', r'socket|tcp_'), ('openmp', r'!\$omp')]
DOWNLOAD = re.compile(r'\bcurl\b|\bwget\b|https?://', re.I)
GENERATED = re.compile(r'simulate_(particles|movie|nanoparticle|noise)|molecule_data|betagal_1jyx|sars_cov2|%simulate|make_random|gauran|ran3|%ran\b|random_number|spiral', re.I)
COMMITTED = re.compile(r"\.txt'|\.pdb'|\.cif'|\.star'|test_data|fixture", re.I)
USER = re.compile(r"(?:checkvar|defined|get_carg)\('(vol1|vol2|stk|projfile|filetab|pdbfile|mskfile|fname|dir_movies|gainref)'|SIMPLE_TEST_VOL1", re.I)
INTRINSIC_CALLS = {'date_and_time', 'cpu_time', 'system_clock', 'random_seed', 'random_number', 'get_command_argument',
                   'get_command', 'get_environment_variable', 'execute_command_line', 'sleep', 'flush', 'move_alloc',
                   'exit', 'abort', 'backtrace', 'srand', 'random_init', 'omp_set_num_threads'}
# launchers that force a separate process; socket and openmp are recorded but do not by themselves decide the tier
ISOLATING = {'coarray', 'mpi', 'offload', 'openacc', 'cuda'}
# router area -> provisional library suite (section 5.2.1 of the plan)
LIB_SUITE = {'fft': 'lib_fft', 'geometry': 'lib_geometry', 'masks': 'lib_masks', 'numerics': 'lib_numerics',
             'optimize': 'lib_optimize', 'stats': 'lib_stats', 'io': 'lib_io', 'class': 'lib_project',
             'utils': 'lib_project', 'parallel': 'lib_parallel', 'network': 'lib_project',
             'highlevel': 'workflow', 'single': 'workflow', 'stream': 'workflow'}
BENCH = re.compile(r'benchmark|\btic\(|\btoc\(|timing|speed|throughput|elapsed', re.I)
WORKFLOW = re.compile(r'%execute\s*\(', re.I)


def read(path):
    with open(path, errors='replace') as fh:
        return fh.read()


def strip_comments(text):
    out = []
    for ln in text.split('\n'):
        q = None
        for i, c in enumerate(ln):
            if q:
                if c == q:
                    q = None
            elif c in '\'"':
                q = c
            elif c == '!':
                ln = ln[:i]
                break
        out.append(ln)
    return '\n'.join(out)


class Test:
    def __init__(self, name):
        self.name = name
        self.routes = {}          # 'standalone' | 'exec' -> dict(files=[...], text=str, lines=int)
        self.area = None
        self.uses = set()
        self.calls = set()        # (object-type-hint, method) or ('', proc)
        self.fail = 'none'
        self.terminators = set()
        self.fixtures = set()
        self.args = set()
        self.launcher = set()
        self.callers = []
        self.timing = {}          # (label, route) -> (status, wall)
        self.twins = []
        self.cluster = []
        self.unique = set()
        self.tier = ''
        self.tier_why = ''


def footprint(text):
    t = strip_comments(text)
    uses = {u.lower() for u in RE_USE.findall(t)} - NOT_PRODUCTION
    uses = {u for u in uses if (u.startswith('simple_') or u.startswith('single_')) and not u.endswith('_tester')}
    # procedures the text defines itself (internal helpers, sub-suite bodies) are not production calls
    local = {m.lower() for m in re.findall(r'^\s*(?:pure\s+|elemental\s+|recursive\s+|logical\s+|integer\s+|real\s+)*(?:subroutine|function)\s+(\w+)', t, re.I | re.M)}
    calls = set()
    for obj, meth in RE_TBP.findall(t) + RE_TBP_BARE.findall(t):
        if obj.lower() in NOT_OBJECT:
            continue
        m = meth.lower()
        if m in ('new', 'kill', 'get', 'set', 'to_char', 'defined', 'get_carg', 'get_rarg', 'get_iarg', 'print', 'write', 'read'):
            continue
        calls.add(('%', m))
    for proc in RE_CALL.findall(t):
        p = proc.lower()
        if p.startswith(('assert_', 'begin_test', 'end_test', 'report_', 'reset_test', 'run_suite', 'simple_end', 'simple_touch',
                         'exec_cmdline', 'del_file', 'simple_mkdir', 'simple_chdir', 'simple_getcwd', 'simple_exception')):
            continue
        if p in INTRINSIC_CALLS or p in local:
            continue
        calls.add(('', p))
    return uses, calls


def classify_failure(text):
    t = strip_comments(text)
    if re.search(r'\bassert_\w+\s*\(|tests_failed|report_summary', t, re.I):
        return 'assertion'
    if re.search(r'\berror\s+stop\b', t, re.I):
        return 'error stop'
    if re.search(r'THROW_HARD', t):
        return 'THROW_HARD'
    return 'none'


def terminators(text):
    t = strip_comments(text)
    out = set()
    if re.search(r'^\s*stop\b', t, re.I | re.M): out.add('stop')
    if re.search(r'\berror\s+stop\b', t, re.I): out.add('error stop')
    if re.search(r'simple_end\s*\(', t, re.I): out.add('simple_end')
    if re.search(r'THROW_HARD', t): out.add('THROW_HARD')
    return out


def fixtures(text):
    t = strip_comments(text)
    f = set()
    if DOWNLOAD.search(t): f.add('download')
    if USER.search(t): f.add('user-supplied')
    if GENERATED.search(t): f.add('generated')
    if COMMITTED.search(t): f.add('committed')
    return f or {'none'}


def launcher(text):
    t = strip_comments(text)
    return {name for name, pat in PLATFORM if re.search(pat, t, re.I)}


# ---------------------------------------------------------------------------
# collection
# ---------------------------------------------------------------------------

def collect_standalone(tests):
    tdir = os.path.join(ROOT, 'production', 'tests')
    for f in sorted(glob.glob(os.path.join(tdir, 'simple_test_*.f90'))):
        name = os.path.basename(f)[len('simple_test_'):-4]
        helpers = sorted(g for g in glob.glob(os.path.join(tdir, 'simple_%s_*.f90' % name)) if g != f)
        text = read(f) + ''.join(read(h) for h in helpers)
        t = tests.setdefault(name, Test(name))
        t.routes['standalone'] = dict(files=[os.path.relpath(x, ROOT) for x in [f] + helpers], text=text, lines=text.count('\n'))
        for k in re.findall(r"checkvar\('(\w+)'", strip_comments(text), re.I):
            t.args.add(k)


def collect_exec(tests):
    # UI: names, areas, required inputs
    required = {}
    for f in glob.glob(os.path.join(ROOT, 'src', 'main', 'ui', 'simple_test', '*.f90')):
        s = read(f)
        for m in re.finditer(r"call\s+(\w+)%new\(\s*&?\s*\n?\s*&?'([^']+)'", s):
            var, name = m.group(1), m.group(2)
            body = s[m.end():]
            nxt = re.search(r"call add_ui_program\(", body)
            body = body[:nxt.start()] if nxt else body
            keys = [a.group(1) for a in re.finditer(r"call\s+%s%%add_input\(\s*UI_\w+,\s*'([^']+)'(.*?)\)\s*\n" % var, body, re.S)
                    if re.search(r",\s*\.true\.\s*,", a.group(2))]
            required[name] = keys
    # commanders: type -> execute procedure, procedure -> body
    bindings, bodies, cfiles = {}, {}, {}
    for f in sorted(glob.glob(os.path.join(ROOT, 'src', 'main', 'commanders', 'test', '*.f90'))):
        s = read(f)
        head = s[:s.lower().find('\ncontains')] if '\ncontains' in s.lower() else ''
        for m in re.finditer(r'type,\s*extends\(commander_base\)\s*::\s*(\w+)(.*?)end\s+type', s, re.S | re.I):
            b = re.search(r'execute\s*=>\s*(\w+)', m.group(2), re.I)
            if b:
                bindings[m.group(1).lower()] = b.group(1).lower()
        for m in RE_SUB.finditer(s):
            bodies[m.group(1).lower()] = (os.path.relpath(f, ROOT), m.group(0), head)
    # routers: case name -> variable -> type
    for f in sorted(glob.glob(os.path.join(ROOT, 'src', 'main', 'exec', 'simple_test_exec_*.f90'))):
        area = os.path.basename(f)[len('simple_test_exec_'):-4]
        s = read(f)
        decl = {v.lower(): t.lower() for t, v in re.findall(r'type\((\w+)\)\s*::\s*(\w+)', s)}
        for m in re.finditer(r"case\s*\(\s*'(\w+)'\s*\)\s*\n\s*call\s+(\w+)%execute", s, re.I):
            name, var = m.group(1), m.group(2).lower()
            proc = bindings.get(decl.get(var, ''), '')
            t = tests.setdefault(name, Test(name))
            t.area = area
            if proc in bodies:
                file, body, head = bodies[proc]
                t.routes['exec'] = dict(files=[file], text=body, lines=body.count('\n'), head=head, proc=proc)
            else:
                t.routes['exec'] = dict(files=[os.path.relpath(f, ROOT)], text='', lines=0, head='', proc=proc or '?')
            for k in required.get(name, []):
                t.args.add(k)


def attach_unit_suites(tests):
    """Fold the sub-suites each unit_<area> commander runs into that identity's footprint.

    The unit_* commanders are a few lines each; the coverage sits in the *_tester
    modules (and local sub-suites) that suites_<area>() registers in
    simple_commanders_test_class.f90. Without this the coverage accounting would
    treat everything merged into a tester module as lost."""
    cls_path = os.path.join(ROOT, 'src', 'main', 'commanders', 'test', 'simple_commanders_test_class.f90')
    if not os.path.exists(cls_path):
        return
    cls = read(cls_path)
    proc2mod = {m.group(2).lower(): m.group(1).lower()
                for m in re.finditer(r'^\s*use\s+(\w+)\s*,\s*only\s*:\s*(\w+)', cls, re.I | re.M)}
    modfiles = {}
    for dp, _, fs in os.walk(os.path.join(ROOT, 'src')):
        for f in fs:
            if f.endswith('_tester.f90'):
                modfiles[f[:-4].lower()] = os.path.join(dp, f)
    for m in re.finditer(r'subroutine\s+suites_(\w+)\s*\(.*?end\s+subroutine\s+suites_\1', cls, re.S | re.I):
        area = m.group(1).lower()
        procs = re.findall(r"call\s+add_suite\(\s*s\s*,\s*n\s*,\s*'[^']*'\s*,\s*(\w+)\s*\)", m.group(0), re.I)
        extra, files = '', []
        for p in procs:
            mod = proc2mod.get(p.lower(), '')
            if mod in modfiles:
                extra += read(modfiles[mod])
                files.append(os.path.relpath(modfiles[mod], ROOT))
            else:
                mm = re.search(r'^\s*subroutine\s+%s\b.*?^\s*end\s+subroutine\s+%s\b' % (p, p), cls, re.S | re.I | re.M)
                if mm:
                    extra += mm.group(0)
        t = tests.get('unit_' + area)
        if t and 'exec' in t.routes:
            t.routes['exec']['text'] += extra
            t.routes['exec']['files'] += files
            t.routes['exec']['suites'] = procs


def collect_callers(tests):
    files = []
    for d in ('.github', 'scripts', 'doc', 'nice'):
        for dp, _, fs in os.walk(os.path.join(ROOT, d)):
            for f in fs:
                if f.endswith(('.yml', '.yaml', '.sh', '.py', '.pl', '.md', '.txt')):
                    files.append(os.path.join(dp, f))
    texts = {f: read(f) for f in files}
    for t in tests.values():
        pat = re.compile(r'simple_test_%s\b|test=%s\b' % (re.escape(t.name), re.escape(t.name)))
        t.callers = sorted(os.path.relpath(f, ROOT) for f, s in texts.items() if pat.search(s))


def analyse(tests):
    for t in tests.values():
        for r, info in t.routes.items():
            text = info['text'] + info.get('head', '')
            u, c = footprint(text)
            t.uses |= u
            t.calls |= c
            t.fixtures |= fixtures(info['text'])
            t.launcher |= launcher(info['text'])
            t.terminators |= terminators(info['text'])
            fp = classify_failure(info['text'])
            order = ['none', 'THROW_HARD', 'error stop', 'assertion']
            if order.index(fp) > order.index(t.fail):
                t.fail = fp
        if 'none' in t.fixtures and len(t.fixtures) > 1:
            t.fixtures.discard('none')
    # unique coverage and overlap
    owners = collections.defaultdict(set)
    for t in tests.values():
        for c in t.calls:
            owners[c].add(t.name)
    for t in tests.values():
        t.unique = {c for c in t.calls if owners[c] == {t.name}}
        if len(t.routes) == 2:
            t.twins.append('same name on both routes')
    names = sorted(tests)
    pairs = []
    for a, b in itertools.combinations(names, 2):
        A, B = tests[a].calls, tests[b].calls
        if len(A) < 3 or len(B) < 3:
            continue
        j = len(A & B) / len(A | B)
        sub = A <= B or B <= A
        if j >= 0.5 or (sub and min(len(A), len(B)) >= 5):
            pairs.append((a, b, j, sub))
    for a, b, j, sub in pairs:
        tests[a].cluster.append('%s (%.0f%%%s)' % (b, 100 * j, ', subset' if sub else ''))
        tests[b].cluster.append('%s (%.0f%%%s)' % (a, 100 * j, ', subset' if sub else ''))
    return pairs


def load_timing(tests, paths):
    for p in paths:
        label = os.path.basename(os.path.dirname(os.path.abspath(p))) or p
        with open(p) as fh:
            next(fh)
            for ln in fh:
                parts = ln.rstrip('\n').split('\t')
                if len(parts) < 5:
                    continue
                route, name, status, wall = parts[0], parts[1], parts[2], parts[3]
                if name in tests:
                    tests[name].timing[(label, route)] = (status, float(wall))


def propose_tier(t):
    """The fast gate is the `units` sub-suites (owner decision 2026-09-22);
    every other identity is proposed for the extensive tier, manual use,
    the platform label, or deletion. A proposal is confirmed or overruled by
    the review verdict."""
    walls = [w for (s, w) in t.timing.values() if s in ('pass', 'fail')]
    wall = max(walls) if walls else None
    text = ''.join(r['text'] for r in t.routes.values())
    if t.name.startswith('unit_'):
        return 'fast', 'area suite of the fast gate (Phase 2)'
    if t.name == 'units':
        return 'fast (umbrella, not registered)', 'runs every area suite in one process; developer convenience'
    if t.name == 'forked_process':
        return 'platform', 'real child processes, clock polling; excluded from the build by decision'
    if t.launcher & ISOLATING:
        return 'platform', 'uses ' + ', '.join(sorted(t.launcher & ISOLATING))
    if 'socket' in t.launcher and re.search(r'socket_(client|server|comm)', t.name):
        return 'platform', 'socket client/server role'
    if 'download' in t.fixtures or 'user-supplied' in t.fixtures:
        return 'manual', 'needs ' + ', '.join(sorted(t.fixtures & {'download', 'user-supplied'}))
    if BENCH.search(text) and t.fail == 'none':
        return 'delete candidate', 'benchmark or timing tool without assertions'
    if WORKFLOW.search(text):
        if t.fail == 'none':
            return 'workflow (needs assertion)', 'runs production commanders but checks nothing'
        return 'workflow', 'runs production commanders'
    lib = LIB_SUITE.get(t.area or '', 'lib_?')
    if t.fail == 'none' and len(t.calls) < 3:
        return 'delete candidate', 'no failure path and fewer than 3 production calls'
    if t.fail == 'none':
        return '%s (needs assertion) or delete' % lib, 'no failure path; %d production calls' % len(t.calls)
    why = 'assertion-bearing' + ('; measured %.1f s' % wall if wall is not None else '')
    return lib, why


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

RUN_STATE = {'pass': 'measured', 'fail': 'measured (failed)', 'timeout': 'timed out', 'crash': 'crashed', 'noargs': 'missing fixture'}


def timing_cell(t):
    """Phase 0 run state per (label, route): a measured time where the test
    could run, otherwise the reason it could not (section 8 of the plan)."""
    if not t.timing:
        if t.launcher & ISOLATING:
            return 'unsupported capability (not run)'
        return 'not run'
    parts = []
    for (label, route), (status, wall) in sorted(t.timing.items()):
        state = RUN_STATE.get(status, status)
        if status in ('pass', 'fail'):
            parts.append('%s/%s %s %.1fs' % (label, route[:4], state, wall))
        elif status in ('crash', 'fail') and (t.launcher & ISOLATING):
            parts.append('%s/%s unsupported capability' % (label, route[:4]))
        else:
            parts.append('%s/%s %s' % (label, route[:4], state))
    return '; '.join(parts)


def dossier(t):
    L = ['# %s' % t.name, '', 'Area: %s  |  proposed tier: **%s** (%s)' % (t.area or 'unassigned', t.tier, t.tier_why), '']
    L.append('## Routes')
    for r, info in sorted(t.routes.items()):
        extra = ('  (procedure `%s`)' % info['proc']) if r == 'exec' else ''
        L.append('- %s: %s, %d lines%s' % (r, ', '.join('`%s`' % f for f in info['files']), info['lines'], extra))
    L += ['', '## Failure path', '', '%s; terminators present: %s' % (t.fail, ', '.join(sorted(t.terminators)) or 'none'), '']
    L += ['## Fixtures, arguments, launcher', '',
          '- fixtures: %s' % ', '.join(sorted(t.fixtures)),
          '- arguments: %s' % (', '.join(sorted(t.args)) or 'none'),
          '- launcher/platform: %s' % (', '.join(sorted(t.launcher)) or 'serial'), '']
    L += ['## Timing', '', timing_cell(t), '']
    L += ['## Footprint', '', 'Production modules imported (%d): %s' % (len(t.uses), ', '.join(sorted(t.uses)) or '-'), '']
    calls = sorted(('%s%s' % ('%' if k == '%' else '', m)) for k, m in t.calls)
    L += ['Production calls (%d): %s' % (len(calls), ', '.join(calls) or '-'), '']
    uniq = sorted(('%s%s' % ('%' if k == '%' else '', m)) for k, m in t.unique)
    L += ['**Unique to this test (%d):** %s' % (len(uniq), ', '.join(uniq) or 'nothing; every call it makes is also made by another test'), '']
    L += ['## Overlap candidates', '']
    L += ['- ' + x for x in t.twins + t.cluster] or ['- none found by footprint']
    L += ['', '## Callers outside the tests', '']
    L += ['- ' + c for c in t.callers] or ['- none (only the test tree names it)']
    L += ['', '## Verdict', '', '_keep | modify | merge into <id> | demote | delete | retire | investigate_ — reviewer, date, note', '']
    return '\n'.join(L)


def inventory(tests, pairs, out):
    by_area = collections.defaultdict(list)
    for t in tests.values():
        by_area[t.area or 'unassigned'].append(t)
    n_s = sum('standalone' in t.routes for t in tests.values())
    n_e = sum('exec' in t.routes for t in tests.values())
    n_both = sum(len(t.routes) == 2 for t in tests.values())
    nofail = sum(t.fail == 'none' for t in tests.values())
    L = ['# Test inventory', '',
         'Generated by `scripts/test_review_dossier.py`; the verdict and note columns are filled in by',
         'reviewers (section 9 of `uniform_test_environment_refactoring.md`) and are the only columns',
         'edited by hand. Re-running the script keeps them.', '',
         '%d canonical test identities from %d route implementations (%d standalone programs, %d `simple_test_exec` cases; %d identities on both routes); %d identities with no failure path on any route; %d footprint-overlap pairs. Terms as defined in section 4 of the plan.'
         % (len(tests), n_s + n_e, n_s, n_e, n_both, nofail, len(pairs)), '',
         'Tier proposals follow section 5 of the plan: `units` is the fast gate; everything else is proposed for a library suite (`lib_<area>`, from its router area), the workflow gates, manual use, the platform label, or deletion, to be confirmed or overruled by the review verdict.', '']
    cols = ['test', 'routes', 'lines', 'failure path', 'run state / time', 'fixtures / args', 'launcher', 'overlap', 'proposed tier', 'verdict', 'note']
    for area in sorted(by_area):
        L += ['## %s' % area, '', '| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
        for t in sorted(by_area[area], key=lambda x: x.name):
            routes = '+'.join(sorted(r[0].upper() for r in t.routes))
            lines = '/'.join(str(t.routes[r]['lines']) for r in sorted(t.routes))
            fx = ', '.join(sorted(t.fixtures)) + ((' / ' + ' '.join(sorted(t.args))) if t.args else '')
            ov = '; '.join(t.twins + t.cluster) or '-'
            row = [t.name, routes, lines, t.fail, timing_cell(t), fx, ', '.join(sorted(t.launcher)) or '-', ov,
                   '%s (%s)' % (t.tier, t.tier_why), t.verdict, t.note]
            L.append('| ' + ' | '.join(str(x).replace('|', '\\|') for x in row) + ' |')
        L.append('')
    L += ['## Retired tests', '', '| test | date | reason | replacement |', '|---|---|---|---|'] + keep_retired(out) + ['']
    with open(out, 'w') as fh:
        fh.write('\n'.join(L))


def keep_retired(path):
    """Rows of the retired-tests table in an existing inventory (hand-written; kept verbatim)."""
    if not os.path.exists(path):
        return []
    rows, active = [], False
    for ln in open(path):
        if ln.startswith('## '):
            active = ln.startswith('## Retired tests')
            continue
        if active and ln.startswith('| ') and not ln.startswith('| test ') and not ln.startswith('|---'):
            rows.append(ln.rstrip('\n'))
    return rows


def keep_verdicts(tests, path):
    """Carry hand-written verdict/note columns over from an existing inventory."""
    for t in tests.values():
        t.verdict, t.note = '', ''
    if not os.path.exists(path):
        return {}
    kept = {}
    for ln in open(path):
        if not ln.startswith('| ') or ln.startswith('| test ') or ln.startswith('|---'):
            continue
        cells = [c.strip() for c in ln.strip().strip('|').split(' | ')]
        if len(cells) >= 11 and cells[0] in tests:
            tests[cells[0]].verdict, tests[cells[0]].note = cells[9], cells[10]
            kept[cells[0]] = cells[9]
    return kept


def coverage_after(tests, path):
    verdicts = keep_verdicts(tests, path)
    dropped = {n for n, v in verdicts.items() if re.match(r'(delete|merge into|retire)\b', v, re.I)}
    before = collections.defaultdict(set)
    after = collections.defaultdict(set)
    for t in tests.values():
        for c in t.calls:
            before[c].add(t.name)
            if t.name not in dropped:
                after[c].add(t.name)
    lost = sorted(c for c in before if not after[c])
    mods_before = set().union(*(t.uses for t in tests.values()))
    mods_after = set().union(*(t.uses for t in tests.values() if t.name not in dropped))
    print('verdicts read: %d; identities dropped by delete/merge/retire: %d' % (len(verdicts), len(dropped)))
    print('production calls exercised before %d, after %d; lost: %d' % (len(before), len(after), len(lost)))
    for c in lost:
        print('   %s%s   (was covered only by: %s)' % ('%' if c[0] == '%' else '', c[1], ', '.join(sorted(before[c]))))
    print('production modules imported by any test before %d, after %d; lost: %s' % (
        len(mods_before), len(mods_after), ', '.join(sorted(mods_before - mods_after)) or 'none'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--timing', nargs='*', default=[], help='results.tsv files from test_timing_run.sh')
    ap.add_argument('--inventory', default=os.path.join(ROOT, 'doc', 'refactoring_notes', 'test_inventory.md'))
    ap.add_argument('--dossiers', default=os.path.join(ROOT, 'build_test_runs', 'dossiers'))
    ap.add_argument('--coverage-after', metavar='INVENTORY', help='report coverage lost by the verdicts in INVENTORY and exit')
    a = ap.parse_args()
    tests = {}
    collect_standalone(tests)
    collect_exec(tests)
    attach_unit_suites(tests)
    collect_callers(tests)
    pairs = analyse(tests)
    load_timing(tests, a.timing)
    for t in tests.values():
        t.tier, t.tier_why = propose_tier(t)
    if a.coverage_after:
        coverage_after(tests, a.coverage_after)
        return 0
    keep_verdicts(tests, a.inventory)
    os.makedirs(a.dossiers, exist_ok=True)
    for t in tests.values():
        with open(os.path.join(a.dossiers, t.name + '.md'), 'w') as fh:
            fh.write(dossier(t))
    inventory(tests, pairs, a.inventory)
    tiers = collections.Counter(t.tier for t in tests.values())
    print('%d identities; %d dossiers in %s; inventory at %s' % (len(tests), len(tests), os.path.relpath(a.dossiers, ROOT), os.path.relpath(a.inventory, ROOT)))
    for k, v in sorted(tiers.items(), key=lambda kv: -kv[1]):
        print('   %-30s %d' % (k, v))
    return 0


if __name__ == '__main__':
    sys.exit(main())
