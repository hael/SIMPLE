#!/usr/bin/env bash

set -euo pipefail

for dir in */; do
  # Skip if not a directory (extra safety)
  [ -d "$dir" ] || continue

  echo "Processing $dir"

  # Move all contents (including hidden files, safely)
  shopt -s dotglob nullglob
  mv "$dir"* .
  shopt -u dotglob nullglob

  # Remove the empty directory
  rmdir "$dir"
done

