#!/usr/bin/env python3
"""check_guide_imports.py -- every @import of a CLAUDE.md under src/ and production/ names a file.

A folder's CLAUDE.md imports its guide (.github/skills/<name>/SKILL.md) so that Claude Code
loads it when a file in the folder is read. An import whose target is missing is dropped
silently, so a renamed or moved skill would stop its guide from loading without notice.
usage: scripts/check_guide_imports.py [ROOT]
exit status: 0 every import resolves, 1 an import is missing
"""

import os
import re
import sys

TARGET_DIRS = ["src", "production"]
IMPORT_RE = re.compile(r"^@(\S+)\s*$")


def main():
    root = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
    missing = []
    for target in TARGET_DIRS:
        for dirpath, dirnames, filenames in os.walk(os.path.join(root, target)):
            if "CLAUDE.md" not in filenames:
                continue
            path = os.path.join(dirpath, "CLAUDE.md")
            with open(path, "r", encoding="utf-8") as f:
                for lineno, line in enumerate(f, 1):
                    m = IMPORT_RE.match(line.strip())
                    if m and not os.path.isfile(os.path.join(dirpath, m.group(1))):
                        missing.append("{}:{}: {}".format(os.path.relpath(path, root), lineno, m.group(1)))
    for entry in missing:
        print("check_guide_imports: missing import target " + entry, file=sys.stderr)
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
