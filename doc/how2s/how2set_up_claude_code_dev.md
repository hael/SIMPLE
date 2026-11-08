# How to set up a machine for SIMPLE development with Claude Code

This recipe sets up a Linux box (or a Mac) where Claude Code runs the full
development cycle on SIMPLE: specify, code, compile, test, iterate. You can
hand it to a Claude Code session on the new machine ("follow
doc/how2s/how2set_up_claude_code_dev.md"). The session does every step it can.
For the steps that need you (sudo, logging in, approving a change to your
shell start-up file), it stops and asks.

For unattended runs of a written plan, see
[how2run_autonomous_agents.md](how2run_autonomous_agents.md) instead. This
recipe covers interactive work, in a terminal or in the VS Code extension.

## How the setup is layered

A Claude Code session in a SIMPLE checkout reads its instructions and settings
from four places. Only the first one comes with `git clone`.

| Layer | Where | What it holds | Who sets it up |
|---|---|---|---|
| Repository | `CLAUDE.md`, the folder `CLAUDE.md` files, `.claude/skills` (a link to `.github/skills`), `.claude/settings.json` | How to work in SIMPLE: skill routing, layering, Fortran and compile-time rules, the build and test cycle, the permissions for build and test commands, `git commit` and `git push` denied | Tracked; nothing to do |
| Build tree | `build/tools/fortls`, `build/claude-plugins`, `.claude/settings.local.json` | The Fortran language server plugin (diagnostics after edits, go-to-definition) | The first build (step 5) |
| Machine | `~/.bashrc`, `~/.claude/CLAUDE.md`, `~/.claude/settings.json` | Toolchain, disks, version-control policy, things the agent must not do, the default model | You and the session (steps 2, 3, 6) |
| Memory | `~/.claude/projects/<checkout path>/memory/` | Facts sessions learn on this machine (timings, quirks, run state) | Sessions write it themselves |

The repository layer is the authority on how to work in SIMPLE. The machine
layer only says what is particular to this box. Don't copy repository rules
into it, because the copies go stale.

## Step 1. Install Claude Code and log in (you)

```bash
curl -fsSL https://claude.ai/install.sh | bash
claude --version
claude doctor
```

The installer puts `claude` in `~/.local/bin`. Make sure that directory is on
`PATH` for every shell (see step 3). Then start `claude` once and run
`/login`. For the other routes (a long-lived token, an API key, Bedrock,
Vertex, Foundry), see Part 1, step 5 of
[how2run_autonomous_agents.md](how2run_autonomous_agents.md). In an
interactive session, `/status` shows which credential is in use.

To work in VS Code, install the Claude Code extension. It uses the same
`~/.claude` and the same checkout settings as the terminal.

## Step 2. Install the build toolchain (you, with sudo)

[doc/installation.md](../installation.md) lists what SIMPLE needs:

- GNU `gcc`, `g++` and `gfortran` 14.2 or newer
- CMake
- an OpenMP build of OpenBLAS (CMake builds one in the build tree when the
  system has none)
- FFTW 3.3 with its double, single and threaded builds
- libcurl and libTIFF
- Python 3.7 or newer for the build scripts and the fast gate

Install `ripgrep` (`rg`) too: the repository instructions search with it.
The agent runs under your account and never uses sudo. It tells you which
packages are missing and you install them.

On enterprise Linux the distribution compiler is usually too old. On Oracle
Linux 8, for example, the system GCC is 8.5, and GCC 15 comes from
`gcc-toolset-15`:

```bash
sudo dnf install -y gcc-toolset-15 gcc-toolset-15-gcc-gfortran gcc-toolset-15-gcc-c++ \
    cmake fftw-devel libcurl-devel libtiff-devel openblas-devel git
sudo dnf install -y ripgrep     # from EPEL on enterprise Linux
```

On macOS, follow the "New Mac Installation" section of
[doc/installation.md](../installation.md) (Homebrew GCC and FFTW; BLAS comes
from Accelerate).

## Step 3. Make the toolchain visible to every shell (the session proposes, you approve)

Claude Code runs commands in a non-interactive bash that it initializes from
your profile. Most `~/.bashrc` files return early for non-interactive shells
(`[[ $- != *i* ]] && return` or a `case $- in` block). Anything the build
needs must therefore come **before** that check. That covers the compiler, the
OpenBLAS override, `TMPDIR`, `~/.local/bin` and SIMPLE's own variables.
Otherwise the agent compiles with the wrong compiler while your own terminal
works.

An example from an Oracle Linux 8 box:

```bash
# at the top of ~/.bashrc, before the interactive-shell check
export PATH="$HOME/.local/bin:$PATH"
# GCC 15 from gcc-toolset-15
[ -f /opt/rh/gcc-toolset-15/enable ] && source /opt/rh/gcc-toolset-15/enable
# OpenBLAS 0.3.15 misdetects recent Intel CPUs; use its AVX2 kernels
export OPENBLAS_CORETYPE=Haswell
# a temp dir on a disk that allows exec and has room
export TMPDIR=$HOME/tmp
# SIMPLE (the build installs into build/)
export SIMPLE_EMAIL="my.name@uni.edu"
export SIMPLE_QSYS="local"
export SIMPLE_PATH=$HOME/src/SIMPLE/build
export PATH=${SIMPLE_PATH}/scripts:${SIMPLE_PATH}/bin:${PATH}
```

Keep only the lines this machine needs. `OPENBLAS_CORETYPE` is for a system
OpenBLAS that misdetects the CPU. `TMPDIR` is for a `/tmp` that is `noexec` or
small. The SIMPLE lines are what `add2.bashrc` adds (see
[doc/installation.md](../installation.md)). Shell changes take effect in a new
Claude Code session.

Check from a session (its Bash tool is the non-interactive shell that
matters):

```bash
gfortran --version      # 14.2 or newer, not the system compiler
cmake --version
rg --version
echo "$TMPDIR"
```

## Step 4. Clone the repository

```bash
mkdir -p ~/src && cd ~/src
git clone https://github.com/hael/SIMPLE.git
cd SIMPLE
ls .claude/skills/ | head        # the link to .github/skills resolves
python3 scripts/check_guide_imports.py .
```

Put the checkout and its build on a disk with room. A Release build tree is
about 1.5 GB, and test and benchmark runs need far more. Start Claude Code
from the checkout root (`cd ~/src/SIMPLE && claude`, or open the folder in VS
Code) so that it loads `CLAUDE.md` and the skills.

## Step 5. The first build (the session)

```bash
./compile_clean.sh          # Release build with tests: configure, make, fast gate, install into build/
```

The script builds the library, the executables and the test programs. It runs
the fast gate (`scripts/run_fast_gate.sh build`) and installs into `build/`
only when the gate passes. Use `./compile_debug.sh` for a debug build. If
`find_package(Python3)` picks an interpreter older than 3.7, pass
`-DPython3_EXECUTABLE=<a newer python3>` to CMake.

The first build also sets up the Fortran language server (option `FORTLS`, on
by default; see `scripts/CMakeLists.txt` and `cmake/InstallFortls.cmake`). It
installs a pinned fortls in `build/tools/fortls` and registers the plugin
`fortls-lsp@simple-build` for this checkout in `.claude/settings.local.json`.
That file is never committed: Claude Code keeps it in your global git ignore
file. Registration needs network access and the `claude` CLI on `PATH` (or
`CLAUDE_CLI` set). Without them the build only warns. **Start a new session
after the first build** to load the plugin. `/plugin` then lists it.

Then check that the installed programs run:

```bash
simple_test_exec test=units      # every fast unit suite in one process
```

## Step 6. The machine's instructions: `~/.claude/CLAUDE.md` (the session drafts, you review)

This file loads into every session on the machine, for every project. It holds
what is true of the box and of how you want agents to work, and nothing the
repository already says. Here is a template. The session fills it in from what
it found in steps 2 to 5:

```markdown
# This machine: <name> (<what it is>)

Machine-specific notes for every project on this box. Repo instructions
(CLAUDE.md, AGENTS.md) still apply; this file adds what is particular to the box.

## Hardware and OS
- <OS and version>, <CPU: model, cores, SIMD (AVX2/AVX-512/NEON)>, <RAM>.
  <GPU or "No GPU: CUDA/GPU builds cannot run here">.
- <Hybrid CPU? e.g. "CPUs 0-7 are P-cores: pin single-thread benchmarks with taskset -c 0-7">

## Toolchain (set up in ~/.bashrc, before its interactive-shell check)
- <compiler and where it comes from; the system compiler not to use>
- <CMake version; the Python CMake should use if the default is too old>
- <BLAS quirks, e.g. OPENBLAS_CORETYPE>

## Filesystem
- <TMPDIR and why>
- <which disks hold sources and builds, which hold datasets and large test runs,
  and which must not be filled>

## Building and testing
- This is a build/test box: run the full cycle (specify, code, compile, test,
  iterate) as the repository's CLAUDE.md describes. Keep build output quiet
  and read only errors and failures.

## Version control
- <your policy, e.g. "Leave every change uncommitted in the working tree, so the
  whole task is one reviewable diff. Never commit, amend, stash, push, create or
  switch branches, rebase, reset or restore. Read-only git is fine. End every
  task by listing the changed and new files (git status).">
- Repos: <paths, remotes; how changes move between machines>

## Things not to do
- Don't use sudo: ask the user to run privileged commands.
- <e.g. "Never reboot or shut down: the disk is encrypted and needs the passphrase at the keyboard.">
```

Two rules for this file:

- Don't restate the repository's rules (build commands, skill routing, test
  tiers). They change with the code, and a copy here goes stale. One
  exception: you may make a repository rule stricter, as the version-control
  policy above does.
- Prefer a short line that says *why* over a long one that says only *what*.
  The agent applies a reason to cases the line does not name.

Optional user-wide settings go in `~/.claude/settings.json`, for example the
default model and effort:

```json
{ "model": "opus", "theme": "auto" }
```

The repository's `.claude/settings.json` already allows the build, gate and
test commands without prompting. Add allowances for this machine only (local
data paths, extra tools) with `/permissions` or the `fewer-permission-prompts`
skill. They go into `.claude/settings.local.json`, not into the tracked file.

## Step 7. Check the cycle end to end (the session)

Make a trivial edit, for example to a comment in a module of `src/main/ori`,
and run the day-to-day loop:

```bash
./recompile.sh                                  # make -C build, fast gate, make install
build/bin/simple_test_exec test=unit_ori        # the touched area's suite
git diff --stat                                 # the change is in the working tree, uncommitted
```

Then revert the edit by hand: the edit is yours to undo, not git's. Record
the timings as a memory, along with anything about this machine that the
next session should not have to rediscover.

## The development cycle on a ready machine

Everything below is in the repository's `CLAUDE.md`, the skills and the
policies. This section only shows how the cycle maps onto them.

1. **Specify.** Write what is to change, and how you will know it works, as a
   note under `doc/implementation_notes/planned/` or
   `doc/refactoring_notes/planned/`. Read the matching policy under
   `doc/policies/` first. For work across subsystems, start with the
   `simple-architecture` skill.
2. **Code.** Choose the narrowest workflow skill (routing in `CLAUDE.md`);
   the folder guides load by themselves as files are read. Follow
   `doc/policies/layering_and_project_ownership_policy.md` and
   `doc/policies/compile_time_policy.md`.
3. **Compile.** Run `./recompile.sh`, which is `make -C build -j<n>`, the fast
   gate, then install. A new source file under `src/` triggers a CMake rerun
   by itself.
4. **Test.** Run the touched area's `simple_test_exec test=unit_<area>`, then
   the high-level entries the change affects (`simple_test_exec test=list`;
   section 1.4 of `doc/policies/test_environment_policy.md`). Start long runs
   in the background.
5. **Iterate.** Repeat steps 2 to 4 until the note's exit criteria hold, and
   update the note with what was built, run and measured.
6. **Review.** The whole task is one diff in the working tree. The
   maintainer reviews it and commits.

## When something goes wrong

- **The agent compiles with the wrong compiler, or `make` cannot find
  `gfortran`.** The toolchain lines are after the interactive-shell check in
  `~/.bashrc`. Move them up and start a new session.
- **No Fortran diagnostics after edits.** The session started before the first
  build registered the plugin. Start a new one. If the build warned that
  `fortls-lsp` was not registered, set `CLAUDE_CLI` to the `claude`
  executable and run `make -C build fortls`.
- **The fast gate fails on `check_guide_imports`.** A folder `CLAUDE.md`
  names a skill that does not exist. Claude Code silently drops a broken
  import, so the gate is what catches it.
- **The gate reports untracked source files.** A new file under `src/` or
  `production/` is not in git yet. That is fine while you work, but list it
  for the maintainer, because a commit without it does not build.
