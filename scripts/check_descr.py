#!/usr/bin/env python3

import os
import sys

TARGET_DIRS = ["src", "production", "scripts"]
EXCLUDED_DIR = os.path.join("src", "extlibs")


def file_has_descr_header(filepath):
    """
    Returns True if the first non-empty line of the file
    starts with '!@descr:', otherwise False.
    """
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                stripped = line.strip()
                if stripped == "":
                    continue
                return stripped.startswith("!@descr:")
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return False

    return False


def file_has_multiline_descr(filepath):
    """
    Returns True if the file carries more than one '!@descr:' line.
    The tag is a one-liner for the code map (like the subject line of a
    commit message); further description goes in plain '!' comments below it.
    """
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            return sum(1 for line in f if line.strip().startswith("!@descr:")) > 1
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return False


def is_excluded(path, base_dir):
    """
    Returns True if path is inside the excluded directory.
    """
    excluded_path = os.path.abspath(os.path.join(base_dir, EXCLUDED_DIR))
    path = os.path.abspath(path)
    return os.path.commonpath([path, excluded_path]) == excluded_path


def find_missing_descr(base_dir):
    missing = []
    multiline = []

    for dirname in TARGET_DIRS:
        search_path = os.path.join(base_dir, dirname)

        if not os.path.isdir(search_path):
            continue

        for root, dirs, files in os.walk(search_path):
            # Remove excluded directory from traversal
            dirs[:] = [
                d for d in dirs
                if not is_excluded(os.path.join(root, d), base_dir)
            ]

            for file in files:
                if file.lower().endswith(".f90"):
                    full_path = os.path.join(root, file)

                    if not is_excluded(full_path, base_dir):
                        if not file_has_descr_header(full_path):
                            missing.append(full_path)
                        elif file_has_multiline_descr(full_path):
                            multiline.append(full_path)

    return missing, multiline


if __name__ == "__main__":
    base_directory = sys.argv[1] if len(sys.argv) > 1 else "."

    results, multiline = find_missing_descr(base_directory)

    if results:
        print("Files missing '!@descr:' header:\n")
        for path in results:
            print(path)
    if multiline:
        print("Files with more than one '!@descr:' line (the tag is a one-liner;"
              " continue in plain '!' comments):\n")
        for path in multiline:
            print(path)
    if results or multiline:
        sys.exit(1)
    else:
        print("All checked .f90 files contain a single '!@descr:' line at the beginning.")
        sys.exit(0)

