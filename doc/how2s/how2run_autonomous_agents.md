# How to run an autonomous Claude Code agent on a Linux box from your Mac

This recipe sets up a Linux box that runs Claude Code unattended on a written
plan, phase by phase, while you control it from your Mac over ssh.

The tool that drives it is `arun`: one shell script and a small engine, kept in
`~/agent_runs` on the Mac, outside any repository. Agent tooling does not
belong in the SIMPLE tree. You take arun from an existing installation (the
home Mac) and make it project-neutral, so that it runs SIMPLE plans and plans
on any other git repository.

What you get:

- You write a plan with phases, exit criteria and a file table that says which
  files each phase may edit. `arun launch` copies your checkout to the box and
  starts the run.
- Claude Code works on a private copy of the checkout. It never commits. After
  each phase the driver writes the phase's raw diff, and at the end the
  complete diff against the commit the copy started from.
- When the agent meets something it cannot decide, it stops and writes a stop
  report. You answer in `rulings.md` and resume.
- Every ten minutes a LaunchAgent on the Mac copies the run's state into
  `~/agent_runs/runs/<run>/live/`. Claude on the Mac can then tell you how
  the run is doing without ssh.

## How the pieces fit

```text
Mac                                            Linux box
~/agent_runs/                                  ~/agent_runs/<run>/
  arun  ---- ssh: launch, resume, stop ---->     engine/   copied by arun at launch/resume
  arun.conf    this Mac's defaults               spec/     run.conf, prompt, phase notes, rulings
  engine/      driver, fence, preflight          repo/     private copy of your checkout
  templates/   starting files per project        state/  logs/  scratch/
  runs/<run>/  run.conf prompt.md phases.md      review/   phase_N.diff, full.diff
               rulings.md                        driver.sh  -> claude -p ...  (one session per phase)
               live/    <-- mirror, every 10 min
               result/  <-- arun fetch
```

The Mac holds everything that defines a run. The box holds only copies, the
run's own state, and the box's Claude credentials.

## What you need

- A Linux box you can leave running and on which you have sudo, on x86-64 or
  ARM64. Claude Code documents Ubuntu 20.04+ and Debian 10+; Oracle Linux 8
  works too.
- A Mac on the same network that reaches the box with ssh.
- A copy of `~/agent_runs` from the home Mac.
- A way to pay for Claude: a subscription, an Anthropic API key, or an
  institutional cloud account. Part 1, step 5 covers each.

The agent runs under your own account on the box, never as root.

## Part 1. Prepare the Linux box

### Step 1. Install the tools

```bash
# RHEL, Rocky, Alma, Oracle Linux
sudo dnf install -y git python3 rsync util-linux coreutils curl
# Ubuntu, Debian
sudo apt-get install -y git python3 rsync util-linux coreutils curl
```

These give bash, git, python3 (3.6 or newer), rsync, flock, timeout, setsid
and curl, which arun needs. Add what the project needs to build and test. For
SIMPLE that is the compiler toolchain, CMake and FFTW described in
`doc/installation.md`.

### Step 2. Keep the box up and the run alive

```bash
sudo systemctl mask sleep.target suspend.target hibernate.target hybrid-sleep.target
sudo loginctl enable-linger "$USER"
```

The first line stops desktop power management from suspending the box in the
middle of a run. The second keeps your processes alive after the last ssh
session closes, on systems where logind ends them at logout. arun already
starts the driver with nohup and setsid, which is enough on most systems.
Lingering is cheap insurance.

A reboot ends the run. Nothing is lost: `arun resume <run>` continues it.

### Step 3. Put the runs on a large local disk

A run keeps its copy of the repository, its builds and its test data under
`~/agent_runs/<run>`. A SIMPLE run can fill tens of GB in a day. If your home
directory is small, or on a network share with a quota, point `~/agent_runs`
to a local disk:

```bash
sudo mkdir -p /scratch/"$USER"/agent_runs && sudo chown "$USER": /scratch/"$USER"/agent_runs
ln -s /scratch/"$USER"/agent_runs ~/agent_runs
```

arun keeps a disk floor (`DISK_MIN_GB`, 20 GB by default): below it the
session may only inspect and clean up.

### Step 4. Install Claude Code

```bash
curl -fsSL https://claude.ai/install.sh | bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
~/.local/bin/claude --version
~/.local/bin/claude doctor
```

The installer puts `claude` in `~/.local/bin` and updates it in the
background. A run can last days, and a version change in the middle of a run
makes its sessions harder to compare. Step 5 turns background updates off.
Update by hand between runs with `claude update`.

### Step 5. Give Claude Code its credentials

Claude Code finds its credentials in its environment. arun's driver reads one
file on the box, `~/.config/arun/env`, and exports what it holds to every
session. arun never copies this file; it stays on the box. Create it first:

```bash
mkdir -p ~/.config/arun
install -m 600 /dev/null ~/.config/arun/env
echo 'DISABLE_AUTOUPDATER=1' >> ~/.config/arun/env
```

Then add one of the routes below. If more than one is set, a cloud provider
wins over an API key, and an API key wins over a subscription. An
`ANTHROPIC_API_KEY` anywhere in your environment on the box silently takes
over from a subscription in headless mode. Keep it out of `~/.bashrc` and
other start-up files.

**A. Subscription, logged in once** (Pro, Max, Team or Enterprise). From the
Mac:

```bash
ssh -t workbox ~/.local/bin/claude
# in Claude Code: /login, open the link on the Mac, paste the code back, then /exit
```

The login is stored in `~/.claude/.credentials.json`. Runs count against the
same usage limits as your own use of Claude. When a limit is hit, the driver
sleeps (`LIMIT_SLEEP`, 30 minutes) and continues.

**B. Subscription, long-lived token.** Run `claude setup-token` on the Mac or
on the box. It prints a long-lived token for headless use and does not store
it. Add it to the file:

```bash
CLAUDE_CODE_OAUTH_TOKEN=<the token>
```

**C. Anthropic API key** (Claude Console). Billing is per token. Create the key
in a workspace of its own with a monthly spend limit, so that no run can
overspend. Set `BUDGET_USD` in a run's `run.conf` to cap each session as well.

```bash
ANTHROPIC_API_KEY=sk-ant-...
```

**D. Amazon Bedrock** (an institutional AWS account):

```bash
CLAUDE_CODE_USE_BEDROCK=1
AWS_PROFILE=claude          # a profile in ~/.aws with Bedrock access
                            # (or AWS_BEARER_TOKEN_BEDROCK=<a Bedrock API key> instead)
AWS_REGION=us-east-1        # optional; otherwise the profile's region
ANTHROPIC_MODEL=<an inference profile ID, such as us.anthropic.claude-...>
```

The AWS account needs access to the Anthropic models enabled in the Bedrock
console, and an IAM policy that allows `bedrock:InvokeModel`,
`bedrock:InvokeModelWithResponseStream`, `bedrock:ListInferenceProfiles` and
`bedrock:GetInferenceProfile`.

**E. Google Vertex AI:**

```bash
CLAUDE_CODE_USE_VERTEX=1
CLOUD_ML_REGION=global                  # or your project's region, such as us-east5
ANTHROPIC_VERTEX_PROJECT_ID=<the Google Cloud project ID>
GOOGLE_APPLICATION_CREDENTIALS=$HOME/.config/gcloud/agent-sa.json
                                        # or run `gcloud auth application-default login` once
```

The project needs the Vertex AI API and the Claude models enabled.

**F. Microsoft Foundry:**

```bash
CLAUDE_CODE_USE_FOUNDRY=1
ANTHROPIC_FOUNDRY_RESOURCE=<the resource name, not a URL>
ANTHROPIC_FOUNDRY_API_KEY=<the key>     # leave out to sign in with Entra ID
ANTHROPIC_MODEL=<your Azure deployment name>
```

Fast mode exists only with a subscription or an Anthropic API key. With a
cloud provider, set `FAST_MODE=0` in `run.conf`. Variable names change between
Claude Code versions. If a route does not answer, check it against the pages
under https://code.claude.com/docs/en/ (authentication, amazon-bedrock,
google-vertex-ai, microsoft-foundry, env-vars).

Check from the Mac that the box answers headless, the way the driver will call
it:

```bash
ssh workbox 'set -a; . ~/.config/arun/env; set +a; ~/.local/bin/claude -p "Reply with OK" --output-format json'
```

The reply is JSON with `"result": "OK"` and `"is_error": false`. In an
interactive session on the box, `/status` shows which credential is in use.

## Part 2. Connect the Mac

Make a key if you have none:

```bash
ssh-keygen -t ed25519 -f ~/.ssh/id_ed25519
```

Give the box a short name in `~/.ssh/config`:

```text
Host workbox
    HostName workbox.example.org        # the box's name or IP address
    User <your user on the box>
    IdentityFile ~/.ssh/id_ed25519
    UseKeychain yes
    AddKeysToAgent yes
    ServerAliveInterval 30
    ServerAliveCountMax 4
```

Install the key, load it into the Mac's agent, and test:

```bash
ssh-copy-id -i ~/.ssh/id_ed25519.pub workbox
ssh-add --apple-use-keychain ~/.ssh/id_ed25519
ssh -o BatchMode=yes workbox true && echo ok
```

The last line must print `ok` without asking for anything. The mirror runs from
a LaunchAgent and cannot type a passphrase.

## Part 3. Install arun on the Mac

### 3.1 Copy it

On the home Mac, pack arun without its runs and its logs:

```bash
tar -C ~ -czf ~/Desktop/arun.tgz --exclude 'agent_runs/runs' --exclude 'agent_runs/var' \
    --exclude 'agent_runs/dell_sync' --exclude 'agent_runs/pose_cont' --exclude '__pycache__' agent_runs
```

Carry `arun.tgz` to the work Mac (AirDrop or a USB stick) and unpack it:

```bash
tar -C ~ -xzf arun.tgz
mkdir -p ~/agent_runs/runs ~/agent_runs/var
echo 'export PATH="$HOME/agent_runs:$PATH"' >> ~/.zshrc
```

All paths in the rest of Part 3 are under `~/agent_runs` on the work Mac.

### 3.2 This Mac's defaults: `arun.conf`

arun falls back to the host `dell` and the checkout `~/src/SIMPLE`. Move these
defaults into a file of their own. Create `~/agent_runs/arun.conf`:

```bash
# This Mac's defaults; a run's run.conf overrides them.
HOST=workbox              # the name from ~/.ssh/config
LOCAL_REPO=~/src/SIMPLE   # the checkout arun copies when run.conf names none
TEMPLATE=simple           # arun new <run> starts from templates/<TEMPLATE>/
```

In `arun`, right after the `SPEC_FILES=(...)` line:

```bash
HOST=; LOCAL_REPO=; TEMPLATE=generic
[ -f "$AR/arun.conf" ] && . "$AR/arun.conf"     # this Mac's defaults
DEFAULT_HOST="$HOST"; DEFAULT_LOCAL_REPO="$LOCAL_REPO"
```

In `load()`, replace `HOST=dell; REMOTE_ROOT="agent_runs/$RUN"; LOCAL_REPO="$HOME/src/SIMPLE"` by

```bash
    HOST="$DEFAULT_HOST"; LOCAL_REPO="$DEFAULT_LOCAL_REPO"; REMOTE_ROOT="agent_runs/$RUN"
```

and add after the line that sources `run.conf`:

```bash
    [ -n "$HOST" ] || die "no HOST in arun.conf or runs/$RUN/run.conf"
```

In `cmd_mirror`, replace `${hosts:-dell}` by `${hosts:-$DEFAULT_HOST}`.

### 3.3 One template per project

```bash
cd ~/agent_runs/templates
mkdir simple generic
mv run.conf prompt.md phases.md rulings.md simple/
cp simple/* generic/
```

Delete the `HOST=` and `LOCAL_REPO=` lines from both `run.conf` templates:
`arun.conf` supplies them, and a run sets them only to differ. Then let
`arun new` take a template. Replace `cmd_new()` in `arun` by:

```bash
cmd_new() {
    local run="${1:?usage: arun new <run> [template]}" tpl="$TEMPLATES/${2:-$TEMPLATE}"
    [ -d "$tpl" ] || die "no template $tpl"
    [ -e "$RUNS/$run" ] && die "runs/$run exists already"
    mkdir -p "$RUNS/$run"
    for f in "${SPEC_FILES[@]}"; do
        sed -e "s#@RUN@#$run#g" "$tpl/$f" > "$RUNS/$run/$f"
    done
    echo "created $RUNS/$run from $tpl: fill in run.conf, prompt.md and phases.md, then: arun launch $run"
}
```

### 3.4 The project's settings in `run.conf`

Several parts of arun assume SIMPLE: the fence reads bare paths in a file table
as paths under `src/main/` and lets a listed `x.f90` also cover
`x_tester.f90`; copying a checkout leaves out `build_*/` trees; the preflight
looks for gfortran and CMake and for other checkouts holding
`compile_debug.sh`; `arun status` lists jobs named `simple_`, `ctest` or
`gfortran`. These become settings. Add to `templates/simple/run.conf`:

```bash
# The project. The base prompt, the fence, the preflight and arun status read these.
PROJECT=SIMPLE                       # named in the base prompt
MAINTAINER="Hans Elmlund"            # who delegated the run and reviews the diff
PATH_BASE=src/main/                  # bare paths in the file table live here (empty: the repository root)
COMPANION=".f90:_tester"             # a listed x.f90 also allows x_tester.f90 (empty: none)
COPY_EXCLUDES="/build/ /build_*/ node_modules/"   # left out when arun copies the checkout
TOOLS="gfortran gcc cmake make"      # recorded by the preflight
CHECKOUT_MARKER=compile_debug.sh     # the preflight lists other checkouts holding this file
JOB_PATTERN="simple_|ctest|make|cmake|gfortran|/usr/bin/time"   # the jobs arun status lists
```

and to `templates/generic/run.conf`:

```bash
# The project. The base prompt, the fence, the preflight and arun status read these.
PROJECT=myproject
MAINTAINER="Your Name"
PATH_BASE=
COMPANION=
COPY_EXCLUDES="/build/ node_modules/ .venv/"
TOOLS="gcc make cmake"
CHECKOUT_MARKER=
JOB_PATTERN="make|cmake|ctest|pytest|python3"
```

In both templates set `FAST_MODE=0` unless the box uses a subscription or an
Anthropic API key.

### 3.5 The fence: `engine/fence.py`

Replace the two constants at the top,

```python
TOP_DIRS = ('src/', 'doc/', 'production/', 'scripts/', 'cmake/', '.github/', '.claude/', 'nice/')
FILE_EXT = re.compile(r'\.(f90|F90|f08|inc|inf|md|txt|py|sh|pl|json|cmake|tsv)$')
```

by one that accepts any file extension:

```python
FILE_EXT = re.compile(r'\.[A-Za-z][A-Za-z0-9_]{0,9}$')
```

Replace `resolve()` by these three functions. The top-level directories now
come from the repository itself, and the base of bare paths from `PATH_BASE`:

```python
_TOP = None


def top_dirs():
    """Top-level directories of the repository; a path starting with one is taken as is."""
    global _TOP
    if _TOP is None:
        try:
            _TOP = tuple(d + '/' for d in git('ls-tree', '-d', '--name-only', 'HEAD').splitlines() if d)
        except Exception:
            _TOP = ()
    return _TOP


def path_base():
    """PATH_BASE of run.conf: where bare paths of the file table live ('' = repository root)."""
    base = conf('PATH_BASE').strip('/')
    return base + '/' if base else ''


def resolve(tok, base_dir):
    """Map a path as the plan's file table writes it to a repository path.

    Paths are under PATH_BASE unless they start with a top-level directory. A bare file
    name sits beside the previous path of the row, or directly under PATH_BASE; whichever
    exists wins, a new file takes the row's directory."""
    tok = tok.strip()
    if tok.startswith(top_dirs()):
        return tok
    if '/' in tok:
        return path_base() + tok
    candidates = []
    if base_dir:
        candidates.append(base_dir.rstrip('/') + '/' + tok)
    candidates.append(path_base() + tok)
    for c in candidates:
        if c in tracked_files() or os.path.exists(os.path.join(repo_root(), c)):
            return c
    return candidates[0]
```

In `rules_for()`, replace the tester rule,

```python
                stem, ext = os.path.splitext(r)
                if ext.lower() == '.f90' and not stem.endswith('_tester'):
                    exact.add(stem + '_tester' + ext)
```

by the `COMPANION` setting:

```python
                comp = conf('COMPANION')         # ".f90:_tester": x.f90 also allows x_tester.f90
                if ':' in comp:
                    cext, csuffix = comp.split(':', 1)
                    stem, ext = os.path.splitext(r)
                    if ext.lower() == cext.lower() and not stem.endswith(csuffix):
                        exact.add(stem + csuffix + ext)
```

Update the docstring at the top of the file to say `PATH_BASE` instead of
`src/main/`.

For SIMPLE file tables the new fence allows the same files as the old one,
plus any backticked names whose extensions the old list did not know (such as
`.csv` or `.pdb`).

### 3.6 The base prompt and the driver

In `engine/prompt_base.md`:

- In the first paragraph, write `the @PROJECT@ repository` for `the SIMPLE
  repository` and `@MAINTAINER@` for `Hans Elmlund`.
- Replace `Then read CLAUDE.md and the skills it routes to for the files this
  phase touches.` by `Then read the repository's agent instructions (CLAUDE.md,
  AGENTS.md or similar) and the skills they route to for the files this phase
  touches.`
- Move the bullet that begins `You are explicitly authorized, overriding the
  CLAUDE.md default, to compile` into `templates/simple/prompt.md`, word for
  word. In the base prompt put in its place: `Builds, tests and program runs on
  this machine are authorized as this run's instructions below say. Before the
  first build of the run, find out how the project is built on this machine
  (compiler, modules, environment of the existing checkouts and their build
  trees) and record it in the journal.`
- Add to `templates/generic/prompt.md` a line that says what a run may build
  and run, for example `You may build with make and run the tests with pytest
  on this machine.`

In `engine/driver.sh`, add the two names to `fill()`:

```bash
        -e "s#@PROJECT@#${PROJECT:-$RUN}#g" -e "s#@MAINTAINER@#${MAINTAINER:-the maintainer}#g"
```

### 3.7 The box's credentials file

In both `engine/driver.sh` and `engine/preflight.sh`, right after the line
that sources `spec/run.conf`:

```bash
# this machine's credentials and provider settings; arun never copies this file
if [ -f "$HOME/.config/arun/env" ]; then set -a; . "$HOME/.config/arun/env"; set +a; fi
```

In `engine/preflight.sh`, record which credential the run will use, after the
line `say "claude: ..."`. The order follows Claude Code's precedence, so the
last match is the one it uses:

```bash
auth="subscription login (~/.claude/.credentials.json)"
[ -n "${CLAUDE_CODE_OAUTH_TOKEN:-}" ] && auth="subscription token (CLAUDE_CODE_OAUTH_TOKEN)"
[ -n "${ANTHROPIC_API_KEY:-}" ]       && auth="Anthropic API key (ANTHROPIC_API_KEY)"
[ -n "${CLAUDE_CODE_USE_FOUNDRY:-}" ] && auth="Microsoft Foundry"
[ -n "${CLAUDE_CODE_USE_VERTEX:-}" ]  && auth="Google Vertex AI"
[ -n "${CLAUDE_CODE_USE_BEDROCK:-}" ] && auth="Amazon Bedrock"
say "auth: $auth"
```

In `engine/driver.sh`, let `limit_hit()` recognize the cloud providers'
throttling as well. Its pattern becomes:

```bash
        grep -qiE 'usage limit|rate.?limit|"api_error_status": ?(429|529)|overloaded|quota|exceeded your|ThrottlingException|RESOURCE_EXHAUSTED|too many requests'
```

A spend cap that is reached also looks like a limit. The driver then sleeps and
retries until `MAX_RESUMES` is used up, which costs nothing. Raise the cap, or
stop the run.

### 3.8 The preflight's tool list

In `engine/preflight.sh`, the `== tools` part becomes:

```bash
say "== tools"
for t in git python3 rsync flock timeout setsid ${TOOLS:-}; do
    printf '%-9s %s\n' "$t" "$(command -v "$t" 2>/dev/null || echo MISSING)" >> "$REPORT"
done
{ for t in ${TOOLS:-}; do "$t" --version 2>&1 | head -n 1; done
  python3 --version 2>&1; git --version; } >> "$REPORT" 2>&1
if type module >/dev/null 2>&1; then
    { echo "-- modules"; module -t avail 2>&1 | head -n 60; } >> "$REPORT" 2>&1
fi
if [ -n "${CHECKOUT_MARKER:-}" ]; then
    { echo "-- other checkouts (read-only for the run)"
      find "$HOME" -maxdepth 4 -name "$CHECKOUT_MARKER" -not -path "$ROOT/*" 2>/dev/null | head -n 20; } >> "$REPORT"
fi
for t in git python3 rsync flock timeout setsid; do
    command -v "$t" >/dev/null || { say "missing required tool: $t"; exit 2; }
done
```

rsync is now required (next step).

### 3.9 The job list and the copy

In `engine/remote.sh`, under `status)`, the job list reads `JOB_PATTERN`:

```bash
    jp="$(. "$ROOT/spec/run.conf" >/dev/null 2>&1; echo "${JOB_PATTERN:-make|cmake|ctest}")"
    ps -u "$USER" -o pid,etime,pcpu,rss,args --sort=-pcpu 2>/dev/null |
        grep -E "$jp" | grep -v -e grep -e claude | cut -c1-170 | head -n 10
```

In `arun`, `copy_to()` and `copy_repo()` fall back to tar where the box has no
rsync. The box has rsync (Part 1, step 1), so delete the tar branches. The
copy in `copy_repo()` becomes:

```bash
    local -a pats=() ex=(); local p
    read -r -a pats <<< "${COPY_EXCLUDES:-/build/ node_modules/}"
    for p in "${pats[@]}"; do ex+=(--exclude "$p"); done
    ssh "${SSH_OPTS[@]}" "$HOST" "mkdir -p '$REMOTE_ROOT/repo'"
    rsync -a -e "ssh ${SSH_OPTS[*]}" "${ex[@]}" --exclude '.DS_Store' \
        --exclude-from="$untracked" "$LOCAL_REPO/" "$HOST:$REMOTE_ROOT/repo/"
```

`read -a` keeps the patterns from being expanded as file names. The
LaunchAgent runs arun with macOS's own `/bin/bash` 3.2, so keep new code to
what bash 3.2 has.

### 3.10 Drop the code for the old run layout

arun still carries code for a run layout older than arun itself: a `kit/`
directory, the `POSE_CONT_ROOT` variable, a tmux session and the `migrate`
command. Runs on this Mac never use it. Delete it from `engine/remote.sh`
(`driver_pid`, `run_pids`, the `layout:` line of `status`, `kit/rulings.md`
in `files`, the tmux lines in `start` and `stop`, the whole `migrate`
command), from `arun` (`remote migrate` in `cmd_resume`, the `kit` branch of
`cmd_rulings`, and in `cmd_mirror` the `legacy` label with the two lines that
remove an older mirror), from `engine/driver.sh` (the `kit` lines of `resume_prompt`),
from `engine/cleanup.py` (the `POSE_CONT_ROOT` mark), and the Notes section of
`README.md`. Replace "the Dell" in the README and in `arun`'s help text by
"the remote machine".

### 3.11 Check and start the mirror

```bash
cd ~/agent_runs
bash -n arun engine/*.sh && python3 -m py_compile engine/*.py && echo syntax ok
arun mirror install
```

`arun mirror install` first checks that ssh to the box works without a prompt.
It then installs the LaunchAgent `local.<your user>.arun.mirror`, which runs
`arun sync` every ten minutes and writes `runs/<run>/live/` for every run that
is not idle.

To let Claude on the Mac report on runs, connect the folder `~/agent_runs` in
the Claude app and ask it to remember where the mirror is: "Runs on workbox
are mirrored in ~/agent_runs/runs/<run>/live/; read synced_at.txt first for
freshness, then status.txt and session_tail.txt." Claude's own sandbox cannot
reach the box, so the mirror is how it sees runs.

## Part 4. A smoke test

Before the first real plan, run a toy plan through the whole chain. It takes
a few minutes and checks the login, the fence, the diffs and the mirror.

Make a tiny repository with a plan:

```bash
mkdir -p ~/src/arun_smoke && cd ~/src/arun_smoke && git init -q
cat > plan.md <<'PLAN'
# Smoke test plan

## Phases

Phase 0. Write `calc.py` with a function `add(a, b)` that returns a + b, and
`test_calc.py`, which checks that add(2, 3) == 5 and prints OK.
Exit: `python3 test_calc.py` prints OK.

Phase 1. Add `mul(a, b)` to `calc.py` and a check of mul(2, 3) == 6 to the test.
Exit: `python3 test_calc.py` prints OK.

### File table

| File | Change | Phases |
|------|--------|--------|
| `calc.py` | new | each |
| `test_calc.py` | new | each |

## Progress

(none yet)
PLAN
git add plan.md && git commit -qm "smoke test plan"
arun new smoke generic
printf '# Phase notes for smoke\n' > ~/agent_runs/runs/smoke/phases.md
```

In `~/agent_runs/runs/smoke/run.conf` set:

```bash
LOCAL_REPO=~/src/arun_smoke
PLAN=plan.md
PLAN_DONE=plan.md
ALLOW_TABLE="### File table"
PHASES="0 1"
PHASE_TIMEOUT=1h
TOOLS=
JOB_PATTERN="python3"
```

Replace `~/agent_runs/runs/smoke/prompt.md` by:

```text
The plan @PLAN@ defines two phases under "Phases", each with its exit check.
You may run python3 on this machine. Record each phase under "Progress" in
the plan.
```

Then run it:

```bash
arun launch smoke        # copies, runs the preflight, starts the driver
arun watch smoke         # follow the session; Ctrl-C stops watching, not the run
arun status smoke        # when it has ended
arun fetch smoke
cd ~/src/arun_smoke && git apply --check ~/agent_runs/runs/smoke/result/review/full.diff && echo applies
```

What to expect. The preflight prints the credential in use, the permission
arguments it chose, `write outside the plan's file table -> REFUSED` and
`git push --dry-run -> REFUSED`, then `preflight passed`. Both phases end
`DONE` within minutes. `full.diff` adds `calc.py`, `test_calc.py`, the
Progress entries and the review report `plan_report.md`.

Afterwards remove the toy run: `ssh workbox rm -rf agent_runs/smoke`, then
`rm -rf ~/agent_runs/runs/smoke ~/src/arun_smoke`.

## Part 5. Daily use

1. Write the plan in the repository: phases with exit criteria, and a file
   table whose first column names the files and whose last column names the
   phases that may edit them (`each` for all).
2. `arun new <run>` (or `arun new <run> generic` for a repository other than
   SIMPLE), then fill in `runs/<run>/run.conf` (at least `PLAN`, `PLAN_DONE`,
   `ALLOW_TABLE`, `PHASES`, and `LOCAL_REPO` if it is not the default),
   `prompt.md` and `phases.md`.
3. Check that the fence reads the file table as you meant it. This prints the
   files phase 0 may edit, with their full repository paths:

   ```bash
   d=$(mktemp -d); mkdir "$d/spec"; cp ~/agent_runs/runs/<run>/run.conf "$d/spec/"
   ln -s ~/src/SIMPLE "$d/repo"          # the run's LOCAL_REPO
   ARUN_ROOT=$d python3 ~/agent_runs/engine/fence.py list 0; rm -rf "$d"
   ```

4. `arun launch <run>`. arun copies the tracked files of `LOCAL_REPO` as they
   are in the working tree, uncommitted edits included.
5. Ask Claude on the Mac how the run is doing, or use `arun status <run>` and
   `arun watch <run>`.
6. On a stop, read the stop report, write your decision in
   `runs/<run>/rulings.md`, then `arun resume <run>`. A ruling for a run that
   is going goes out with `arun rulings <run>`; the session reads it when the
   next phase starts.
7. When all phases are done, `arun fetch <run>` and review
   `runs/<run>/result/review/full.diff`. It applies to the commit in
   `result/base_head`.

`arun stop <run>` ends the driver, its Claude session and every job they
started. `arun list` shows all runs from the mirror, without ssh.

## When something goes wrong

- **The preflight says Claude Code does not answer.** The credential is
  missing or wrong. Repeat the check at the end of Part 1, step 5. For a
  subscription login, log in again with `ssh -t workbox ~/.local/bin/claude`.
- **The preflight says the fence hook did not block a write.** python3 is
  missing on the box, or the engine copy is incomplete. Check
  `logs/preflight.txt` in the run directory on the box.
- **The mirror goes stale** (`live/synced_at.txt` says FAILED, or is old). The
  reason is in `live/sync_error.txt`. Most often the ssh key is not in the
  Mac's agent after a restart: `ssh-add --apple-use-keychain`. A VPN on the Mac
  or a network change can also cut the route ("Operation timed out"). The run
  goes on regardless; only the mirror stops.
- **The run waits on a usage limit.** The driver sleeps 30 minutes and tries
  again, up to `MAX_RESUMES` times. Nothing to do unless it is a spend cap.
- **The disk fills.** Below `DISK_MIN_GB` the session may only clean up.
  `arun status <run>` shows what holds the space, inside and outside the run.
- **The box rebooted.** `arun resume <run>` continues the phase in its own
  session.
- **A phase stopped with "file fence or index check failed".** The phase
  changed a file its row in the file table does not allow, perhaps through the
  shell. `logs/phase_N_fence.txt` names the files. Add the row with a ruling,
  or tell the session to undo the change, then resume.

## What the guard rails do, and what they do not

The fence (a hook that runs before every edit and shell command) refuses
edits outside the current phase's files and the run's own directories, any git
command that writes, sudo, ssh, pkill, crontab and systemctl, and removing
`/`, the home directory or `.git`. After each phase the driver checks the
changed files against the file table again, from the diff, so an edit made
through the shell is caught too. The review diff holds everything that differs
from the starting commit, so nothing hides from review. The disk guard keeps
the file system from filling.

These are rails against mistakes, not a sandbox. The session runs as you: it
can read anything you can read, including the credentials file, and the box
keeps its network. The hook does not watch reads, and it matches shell
commands by pattern, which a determined process could get around. Two habits
keep the risk small. Keep a password on sudo (the default), so that the
fence's refusal of sudo is not the only thing between a session and root. And
give the box a credential with its own spend limit.
