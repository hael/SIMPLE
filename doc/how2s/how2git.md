# Useful Git Commands

## Add New Branch To Remote Repo
First, create a new branch locally
```shell
git checkout -b elbrus
```
Then, push it
```shell
git push git@github.com:hael/SIMPLE.git elbrus
```
Remove the local branch with
```shell
git checkout master
git branch -d elbrus
```

## How to Override Local Changes
```shell
git fetch –all git reset –hard origin/master
```

## Change the Editor for Commit Messages
```shell
git config –global core.editor “vi”
```

## Clear All Stashed Branches
```shell
git stash clear
```

## Delete All History of Big Files
```shell
git filter-branch --force --index-filter 'git rm --cached --ignore-unmatch web/SIMPLE2.1/1.0/binaries/simple_linux_120521.tar.gz' --prune-empty --tag-name-filter cat -- --all
```

## Permission Reset
```shell
git config --global --add alias.permission-reset '!git diff -p -R --no-ext-diff --no-color | grep -E "^(diff|(old|new) mode)" --color=never | git apply'
```
