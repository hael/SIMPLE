#!/usr/bin/env python3
# Check the dimensions of all movie files
import sys
import mrcfile, collections
from pathlib import Path

if len(sys.argv) < 2:
    raise SystemExit(f"Usage: {sys.argv[0]} <filetab.txt>")

filetab = sys.argv[1]
dims = collections.defaultdict(list)  # key -> [(entry_idx, filename)]

for idx, line in enumerate(Path(filetab).read_text().splitlines(), start=1):
    f = line.strip()
    if not f:
        continue
    with mrcfile.open(f, permissive=True) as m:
        shp = m.data.shape
    #print(f"{idx}: {f} {shp}")
    if len(shp) == 2:
        key = (shp[1], shp[0], 1)
    else:
        key = (shp[2], shp[1], shp[0])
    dims[key].append((idx, f))

print("\n=== Dimension groups ===")
for k, v in dims.items():
    print(k, len(v))

if not dims:
    raise SystemExit("No readable entries found.")

major_key, major_list = max(dims.items(), key=lambda kv: len(kv[1]))
print(f"\nMajority dimension: {major_key} (count={len(major_list)})")

diff = [(k, items) for k, items in dims.items() if k != major_key]
if not diff:
    print("No different files: all entries match the majority.")
else:
    print("\n=== Different from majority ===")
    for k, items in diff:
        print(f"Dimension {k} (count={len(items)}):")
        for idx, f in items:
            print(f"  entry #{idx}: {f}")
