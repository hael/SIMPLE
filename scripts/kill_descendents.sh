#!/bin/bash
# -*- mode: shell-script -*-
##
## Terminate descendents of the argument process
## SIMPLE runs in batch mode with several parts, in some cases a child
##
## Michael Eager (michael.eager@monash.edu)

E_BADARGS=85

if [ ! -n "$1" ];
then
  echo "Usage: `basename $0` argument1 argument2 etc."
  exit $E_BADARGS
fi


kill_descendant_processes() {
    local pid="$1"
    local and_self="${2:-false}"
    if children="$(pgrep -P "$pid")"; then
        for child in $children; do
            kill_descendant_processes "$child" true
        done
    fi
    if [[ "$and_self" == true ]]; then
        kill -9 "$pid"
    fi
}

if [ $# -ge 1 ]; then
    kill_descendant_processes $1
fi
