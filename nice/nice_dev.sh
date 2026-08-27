#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)
repo_dir=$(cd -- "$script_dir/.." >/dev/null 2>&1 && pwd)
build_dir=${SIMPLE_BUILD_DIR:-"$repo_dir/build"}
python_bin="$build_dir/nice/venv/bin/python3"
manage_py="$script_dir/manage.py"
listen_addr=${NICE_DEV_ADDR:-"127.0.0.1:8000"}

show_help() {
    printf '%s\n' \
        "Usage: nice/nice_dev.sh [runserver options]" \
        "" \
        "Launch NICE directly from the source tree using the existing build" \
        "virtual environment. SIMPLE does not need to be rebuilt after NICE edits." \
        "" \
        "Environment variables:" \
        "  SIMPLE_BUILD_DIR  SIMPLE build directory (default: <repository>/build)" \
        "  NICE_DEV_ADDR     Django listen address (default: 127.0.0.1:8000)" \
        "" \
        "Example:" \
        "  NICE_DEV_ADDR=0.0.0.0:8001 nice/nice_dev.sh"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

if [[ ! -x "$python_bin" ]]; then
    printf 'NICE Python environment not found: %s\n' "$python_bin" >&2
    printf 'Build SIMPLE once with NICE enabled, or set SIMPLE_BUILD_DIR.\n' >&2
    exit 1
fi

if [[ ! -f "$manage_py" ]]; then
    printf 'NICE manage.py not found: %s\n' "$manage_py" >&2
    exit 1
fi

export DJANGO_SETTINGS_MODULE=nice_project.settings_local
export SIMPLE_PATH="$build_dir"
export PATH="$build_dir/bin:$PATH"
export LD_LIBRARY_PATH="$build_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd "$script_dir"

printf 'Preparing NICE development database\n'
"$python_bin" "$manage_py" migrate --noinput
"$python_bin" "$manage_py" ensurelocaluser
"$python_bin" "$manage_py" resetsubmissiontemplate "$build_dir"

printf 'Starting NICE from %s\n' "$script_dir/nice_lite"
printf 'URL: http://%s\n' "$listen_addr"
printf 'Username: %s\n' "${USER:-${LOGNAME:-unknown}}"
printf 'Password: simple\n'

exec "$python_bin" "$manage_py" runserver "$listen_addr" "$@"
