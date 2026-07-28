# Canonical Origin Push Workflow

This applies when `origin` is the canonical `hael/SIMPLE` repository. SIMPLE
practices extreme programming, so changes should be focused, reviewed, and
validated before committing. Do not automatically push commits to
`origin/master`; push only when the user explicitly requests it.

## Before editing

```bash
git branch --show-current
git status --short --branch
git fetch origin
```

Preserve unrelated dirty-worktree changes. Do not reset, clean, or overwrite
them. Read repository instructions and applicable workflow notes first.

## Before committing

```bash
git diff --check
git diff --stat
git diff
```

Run the narrowest relevant validation and report tests that were not run. Do
not claim Linux, BOX, or data tests passed without observing their output.

## Commit

Stage only intended files, inspect the staged diff, then commit:

```bash
git add <intended-files>
git diff --cached --check
git diff --cached --stat
git commit -m "<focused message>"
git status --short --branch
```

## Push only when requested

After explicit user authorization, synchronize with the remote before pushing.
Do not rebase with unrelated uncommitted changes in the worktree:

```bash
git status --short
git fetch origin
git rebase origin/master
```

If the rebase reports conflicts, stop and resolve them deliberately. After
resolving conflicts, run the relevant validation again and confirm the staged
and committed diff before continuing. Do not use force-push to bypass a
rejected push.

Once the rebase succeeds, verify the remote and branch again, then push:

```bash
git remote get-url origin
git branch --show-current
git status --short --branch
git push origin master
```

Stop if the branch is not `master`, the remote is not canonical, or unexpected
files are staged. Never force-push. If `origin/master` advanced and the push is
rejected, fetch, review the divergence, and resolve it explicitly.

## Remote identity check

```bash
git remote get-url origin
git branch --show-current
```

For this workflow, `origin` should be `https://github.com/hael/SIMPLE.git` and
the current branch should be `master`. For fork-based work, use
`.codex/fork-workflow.md` instead.
