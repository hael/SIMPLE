#!/usr/bin/env python3
"""The class-average CNN from start to end with one script: environment, simulated training set from EMDB maps,
training with cross-validation, prediction on your class averages, and a summary report.

  python3 sscnn.py setup   [--torch auto|cuda|cpu]
  python3 sscnn.py build   --set SET_DIR --simple SIMPLE_BUILD [--jobs 4] [--nthr 14] [--entries entries.tsv]
  python3 sscnn.py train   --set SET_DIR --model-dir MODEL_DIR [--recollect]
  python3 sscnn.py predict --model-dir MODEL_DIR --stk cavgs_iterNNN.mrc [--cls classdoc_ranked.txt] [--name NAME]
  python3 sscnn.py report  --model-dir MODEL_DIR [--set SET_DIR]
  python3 sscnn.py all     --set SET_DIR --simple SIMPLE_BUILD --model-dir MODEL_DIR [--stk cavgs.mrc [--cls doc] ...]

Part of SIMPLE (scripts/sscnn/; installed with SIMPLE into $SIMPLE_PATH/scripts/sscnn/). Guide:
doc/how2s/how2train_cavg_cnn.md; the model: CNN.md beside this script.

Steps:
  setup    a virtual environment with numpy and PyTorch (--venv, or $SSCNN_VENV, default
           ~/.cache/simple/sscnn_venv): the CUDA build when nvidia-smi finds a GPU, else the CPU build (--torch to
           choose); skipped when it already imports both. Every other step runs it first.
  build    the simulated set (build_emdb_set.py): per entry of --entries, both half maps and the model downloaded from
           EMDB and RCSB, the map prepared, low-passed and phase-randomised variants at 15, 30 and 50 A, SIMPLE's
           particle simulation, junk, solve2D and the ground truth; restartable (finished runs are kept). Needs SIMPLE
           (--simple, else $SIMPLE_PATH, else the SIMPLE install this script sits in, else the repository's build/:
           the folder holding bin/simple_exec) and internet access.
  train    the classes collected into MODEL_DIR/data (cnn_data.py, once; --recollect to collect again), then
           cross-validation by protein and the final model (cnn_ss.py train) into MODEL_DIR
  predict  the model on class averages (each --stk with its _even and _odd halves; --cls with SIMPLE's class document
           adds populations and resolutions) into MODEL_DIR/predictions/NAME.{tsv,html,json}
  report   MODEL_DIR/summary.html: the set, the training data, the model, the cross-validation and every prediction,
           linked to the detailed pages (cv_report.html, the prediction pages)
  all      setup, build, train, predict (when --stk is given) and report in turn

Options after -- go to the underlying script (build_emdb_set.py, cnn_ss.py train or cnn_ss.py predict), e.g.
  python3 sscnn.py train --set S --model-dir M -- --epochs 30
With `all`, use --build-args, --train-args and --predict-args (each one quoted string).
Each step's output also goes to a log: SET_DIR/build.log, MODEL_DIR/cnn.log. For a long run on a remote machine:
  nohup python3 sscnn.py all ... > sscnn.log 2>&1 &
Standard library only: runs with any Python 3.9 or newer, which then also makes the environment.
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TORCH_CPU_INDEX = 'https://download.pytorch.org/whl/cpu'
DEFAULT_VENV = os.environ.get('SSCNN_VENV') or os.path.join(
    os.environ.get('XDG_CACHE_HOME') or os.path.join(os.path.expanduser('~'), '.cache'), 'simple', 'sscnn_venv')


def say(msg):
    print(f'[sscnn {time.strftime("%H:%M:%S")}] {msg}')
    sys.stdout.flush()


def run_logged(cmd, log, env=None):
    """cmd with its output on the screen and appended to log; its exit code."""
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, 'a') as f:
        f.write(f'\n# {time.strftime("%Y-%m-%d %H:%M:%S")} {" ".join(shlex.quote(c) for c in cmd)}\n')
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, bufsize=1)
        for line in p.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            f.write(line)
            f.flush()
        return p.wait()


def check(rc, what):
    if rc != 0:
        raise SystemExit(f'[sscnn] {what} failed (exit {rc}); see its log')


# ------------------------------------------------------------------ steps
def venv_python(a):
    return os.path.join(os.path.abspath(a.venv), 'bin', 'python')


def setup(a):
    if sys.version_info < (3, 9):
        raise SystemExit('Python 3.9 or newer is needed')
    py = venv_python(a)
    if not os.path.exists(py):
        say(f'creating the environment {a.venv} with {a.python or sys.executable}')
        subprocess.check_call([a.python or sys.executable, '-m', 'venv', a.venv])
    ok = subprocess.run([py, '-c', 'import numpy, torch'], capture_output=True).returncode == 0
    if not ok or a.reinstall:
        want = a.torch
        if want == 'auto':
            want = 'cuda' if shutil.which('nvidia-smi') and subprocess.run(['nvidia-smi', '-L'], capture_output=True).returncode == 0 \
                else 'cpu'
        if sys.platform == 'darwin':
            want = 'default'                                  # macOS wheels: CPU and MPS
        say(f'installing numpy and PyTorch ({want}) into {a.venv}')
        pip = [py, '-m', 'pip', 'install', '--quiet', '--no-cache-dir']
        subprocess.check_call(pip + ['--upgrade', 'pip'])
        subprocess.check_call(pip + ['numpy'])
        subprocess.check_call(pip + ['torch'] + (['--index-url', TORCH_CPU_INDEX] if want == 'cpu' else []))
    r = subprocess.run([py, '-c', 'import sys, numpy, torch; d = ("cuda: " + torch.cuda.get_device_name(0)) if torch.cuda.is_available() '
                        'else ("mps" if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available() else "cpu"); '
                        'print(f"python {sys.version.split()[0]}, numpy {numpy.__version__}, torch {torch.__version__}, device {d}")'],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f'[sscnn] the environment does not import numpy and torch:\n{r.stderr[-2000:]}')
    say(f'environment {a.venv}: {r.stdout.strip()}')
    return py


def find_simple(a):
    """--simple, else $SIMPLE_PATH, else the SIMPLE install holding this script (<prefix>/scripts/sscnn), else the
    repository's build/ (<repo>/scripts/sscnn)."""
    for sp in (a.simple, os.environ.get('SIMPLE_PATH'), os.path.join(HERE, '..', '..'), os.path.join(HERE, '..', '..', 'build')):
        if sp and os.access(os.path.join(sp, 'bin', 'simple_exec'), os.X_OK):
            return os.path.abspath(sp)
        if sp and sp in (a.simple, os.environ.get('SIMPLE_PATH')):
            raise SystemExit(f'[sscnn] no {os.path.join(sp, "bin", "simple_exec")}')
    raise SystemExit('[sscnn] SIMPLE not found: give --simple, the SIMPLE build or install folder holding bin/simple_exec')


def simple_env(a):
    sp = find_simple(a)
    exe = os.path.join(sp, 'bin', 'simple_exec')
    if not os.access(exe, os.X_OK):
        raise SystemExit(f'[sscnn] no {exe}')
    env = dict(os.environ, SIMPLE_PATH=sp)
    env.setdefault('SIMPLE_QSYS', 'local')
    env.setdefault('SIMPLE_EMAIL', 'nobody@example.org')
    r = subprocess.run(['git', '-C', os.path.dirname(sp), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True)
    say(f'SIMPLE {sp}, commit {r.stdout.strip() if r.returncode == 0 else "unknown (not a git checkout)"}')
    return env


def build(a, py, extra):
    env = simple_env(a)
    s = os.path.abspath(a.set)
    os.makedirs(s, exist_ok=True)
    cmd = [py, '-I', os.path.join(HERE, 'build_emdb_set.py'), '--entries', os.path.abspath(a.entries), '--set', s,
           '--jobs', str(a.jobs), '--nthr', str(a.nthr)] + extra
    say(f'building the set in {s} ({a.jobs} entries at a time, {a.nthr} threads each)')
    check(run_logged(cmd, os.path.join(s, 'build.log'), env), 'build')


def train(a, py, extra):
    md = os.path.abspath(a.model_dir)
    data = os.path.join(md, 'data')
    log = os.path.join(md, 'cnn.log')
    if a.recollect or not os.path.exists(os.path.join(data, 'info.json')):
        say(f'collecting the classes of {a.set} into {data}')
        shutil.rmtree(data, ignore_errors=True)
        check(run_logged([py, '-I', os.path.join(HERE, 'cnn_data.py'), '--set', os.path.abspath(a.set), '--out', data], log),
              'collecting the classes')
    say(f'training into {md}')
    check(run_logged([py, '-I', os.path.join(HERE, 'cnn_ss.py'), 'train', '--data', data, '--out', md] + extra, log), 'training')


def predict(a, py, extra):
    md = os.path.abspath(a.model_dir)
    model = os.path.abspath(a.model or os.path.join(md, 'model.pt'))
    if not os.path.exists(model):
        raise SystemExit(f'[sscnn] no model {model}: train first')
    out = os.path.join(md, 'predictions')
    os.makedirs(out, exist_ok=True)
    stks = a.stk or []
    clss = a.cls or []
    names = a.name or []
    if clss and len(clss) != len(stks):
        raise SystemExit('[sscnn] give one --cls per --stk, or none')
    for k, stk in enumerate(stks):
        stk = os.path.abspath(stk)
        name = names[k] if k < len(names) else \
            f'{os.path.basename(os.path.dirname(os.path.dirname(stk)))}_{os.path.basename(os.path.dirname(stk))}_{os.path.splitext(os.path.basename(stk))[0]}'
        cmd = [py, '-I', os.path.join(HERE, 'cnn_ss.py'), 'predict', '--model', model, '--stk', stk,
               '--out', os.path.join(out, f'{name}.tsv'), '--html', os.path.join(out, f'{name}.html')] + \
            (['--cls', os.path.abspath(clss[k])] if clss else []) + extra
        say(f'predicting {stk}')
        check(run_logged(cmd, os.path.join(md, 'cnn.log')), f'prediction of {stk}')


def report(a, py):
    md = os.path.abspath(a.model_dir)
    cmd = [py, '-I', os.path.join(HERE, 'cnn_ss.py'), 'summary', '--model-dir', md] + \
        (['--set', os.path.abspath(a.set)] if a.set else [])
    check(run_logged(cmd, os.path.join(md, 'cnn.log')), 'report')


# ------------------------------------------------------------------ command line
def main():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--venv', default=DEFAULT_VENV, help='the environment {$SSCNN_VENV or ~/.cache/simple/sscnn_venv}')
    common.add_argument('--python', help='the interpreter for a new environment {this one}')
    common.add_argument('--torch', choices=('auto', 'cuda', 'cpu'), default='auto', help='PyTorch build for a new environment {auto}')
    common.add_argument('--reinstall', action='store_true', help='install numpy and PyTorch again')
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='step', required=True)

    def add(name, *groups):
        p = sub.add_parser(name, parents=[common], formatter_class=argparse.RawDescriptionHelpFormatter)
        if 'set' in groups:
            p.add_argument('--set', required=name in ('build', 'train', 'all'), help='the simulated set\'s folder')
        if 'build' in groups:
            p.add_argument('--simple', help='SIMPLE build or install folder holding bin/simple_exec '
                                            '{$SIMPLE_PATH, else the install holding this script, else the repo\'s build/}')
            p.add_argument('--entries', default=os.path.join(HERE, 'entries.tsv'), help='EMDB entries {entries.tsv}')
            p.add_argument('--jobs', type=int, default=4, help='entries at once {4}')
            p.add_argument('--nthr', type=int, default=14, help='threads per SIMPLE program {14}')
        if 'model' in groups:
            p.add_argument('--model-dir', required=True, help='training output: model, reports, predictions')
        if 'train' in groups:
            p.add_argument('--recollect', action='store_true', help='collect the classes again (after the set grew)')
        if 'predict' in groups:
            p.add_argument('--stk', action='append', required=name == 'predict', help='class averages (MRC); repeatable')
            p.add_argument('--cls', action='append', help="SIMPLE's class document for each --stk (classdoc_ranked.txt)")
            p.add_argument('--name', action='append', help='name of each prediction {from the stack\'s path}')
            p.add_argument('--model', help='a model file {MODEL_DIR/model.pt}')
        return p
    add('setup')
    add('build', 'set', 'build')
    add('train', 'set', 'model', 'train')
    add('predict', 'model', 'predict')
    add('report', 'set', 'model')
    al = add('all', 'set', 'build', 'model', 'train', 'predict')
    for k in ('build', 'train', 'predict'):
        al.add_argument(f'--{k}-args', default='', help=f'options for the {k} step (one quoted string)')
    argv = sys.argv[1:]
    extra = argv[argv.index('--') + 1:] if '--' in argv else []
    a = ap.parse_args(argv[:argv.index('--')] if '--' in argv else argv)
    t0 = time.time()
    py = setup(a)
    if a.step == 'build':
        build(a, py, extra)
    elif a.step == 'train':
        train(a, py, extra)
    elif a.step == 'predict':
        predict(a, py, extra)
    elif a.step == 'report':
        report(a, py)
    elif a.step == 'all':
        build(a, py, shlex.split(a.build_args))
        train(a, py, shlex.split(a.train_args))
        if a.stk:
            predict(a, py, shlex.split(a.predict_args))
        report(a, py)
    say(f'{a.step} done in {(time.time() - t0) / 60:.1f} min')


if __name__ == '__main__':
    main()
