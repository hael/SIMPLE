#!/usr/bin/env python3
"""Report per-2D-class defocus statistics from a SIMPLE project file.

This script queries particle metadata via SIMPLE's print_project_vals program
and computes mean/variance of defocus values for particles in each 2D class.

By default:
- only active particles are included (state > 0)
- defocus value per particle is (dfx + dfy) / 2
- population variance is reported

Examples:
  python scripts/defocus_stats_by_cls2d.py my_project.simple
  python scripts/defocus_stats_by_cls2d.py my_project.simple --defocus-field dfx
  python scripts/defocus_stats_by_cls2d.py my_project.simple --include-rejected --sample-variance
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Iterable, Optional

try:
    import matplotlib.pyplot as plt
except ImportError:  # Optional dependency; only required when --plot is used.
    plt = None


def run_cmd(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=False)


def parse_int_token(token: str) -> Optional[int]:
    try:
        return int(token)
    except ValueError:
        try:
            return int(float(token))
        except ValueError:
            return None


def parse_float_token(token: str) -> Optional[float]:
    try:
        return float(token)
    except ValueError:
        return None


def parse_rows(stdout: str, expected_payload_cols: int) -> Iterable[tuple[int, list[str]]]:
    """Yield (state, payload_columns) from print_project_vals output.

    SIMPLE output rows are expected to begin with:
      <index> <state> <payload...>
    Non-data/log lines are ignored.
    """
    for raw in stdout.splitlines():
        line = raw.strip()
        if not line:
            continue
        cols = line.split()
        if len(cols) < 2 + expected_payload_cols:
            continue
        idx = parse_int_token(cols[0])
        state = parse_int_token(cols[1])
        if idx is None or state is None:
            continue
        payload = cols[2 : 2 + expected_payload_cols]
        yield state, payload


def compute_stats(values: list[float], sample_variance: bool) -> tuple[float, float]:
    n = len(values)
    if n == 0:
        return math.nan, math.nan
    mean = sum(values) / n
    if n == 1:
        return mean, 0.0
    ss = 0.0
    for v in values:
        d = v - mean
        ss += d * d
    denom = (n - 1) if sample_variance and n > 1 else n
    return mean, ss / denom


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compute mean/variance of particle defocus values within each 2D class "
            "from a SIMPLE project file."
        )
    )
    parser.add_argument("projfile", type=Path, help="Path to input .simple project file")
    parser.add_argument(
        "--simple-private-exec",
        default="simple_private_exec",
        help="Path to SIMPLE private executable that provides print_project_vals",
    )
    parser.add_argument(
        "--defocus-field",
        choices=["avg", "dfx", "dfy"],
        default="avg",
        help="Defocus value used per particle: avg=(dfx+dfy)/2, or dfx, or dfy",
    )
    parser.add_argument(
        "--include-rejected",
        action="store_true",
        help="Include particles with state <= 0 (default: excluded)",
    )
    parser.add_argument(
        "--sample-variance",
        action="store_true",
        help="Use sample variance (N-1 denominator) instead of population variance",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Create a scatter plot of mean defocus (x) vs variance (y)",
    )
    parser.add_argument(
        "--plot-out",
        type=Path,
        default=Path("defocus_mean_vs_variance.png"),
        help="Output image path used when --plot is enabled",
    )
    args = parser.parse_args()

    if not args.projfile.exists():
        print(f"error: project file not found: {args.projfile}", file=sys.stderr)
        return 2

    keys = ["class", "dfx", "dfy"]
    cmd = [
        args.simple_private_exec,
        "prg=print_project_vals",
        f"projfile={args.projfile}",
        "oritype=ptcl2D",
        f"keys={','.join(keys)}",
    ]
    out = run_cmd(cmd)
    if out.returncode != 0:
        print("error: print_project_vals failed", file=sys.stderr)
        if out.stderr.strip():
            print(out.stderr.strip(), file=sys.stderr)
        return out.returncode or 1

    per_class: DefaultDict[int, list[float]] = defaultdict(list)

    for state, payload in parse_rows(out.stdout, expected_payload_cols=3):
        if not args.include_rejected and state <= 0:
            continue

        cls_id = parse_int_token(payload[0])
        dfx = parse_float_token(payload[1])
        dfy = parse_float_token(payload[2])
        if cls_id is None or dfx is None or dfy is None:
            continue
        if cls_id <= 0:
            continue

        if args.defocus_field == "dfx":
            defocus = dfx
        elif args.defocus_field == "dfy":
            defocus = dfy
        else:
            defocus = 0.5 * (dfx + dfy)

        per_class[cls_id].append(defocus)

    if not per_class:
        print("No matching ptcl2D rows found (after filtering).", file=sys.stderr)
        return 1

    print("class  n_particles  mean_defocus  variance_defocus")
    rows: list[tuple[int, int, float, float]] = []
    for cls_id in sorted(per_class):
        vals = per_class[cls_id]
        mean, var = compute_stats(vals, sample_variance=args.sample_variance)
        rows.append((cls_id, len(vals), mean, var))
        print(f"{cls_id:5d}  {len(vals):11d}  {mean:12.6f}  {var:16.6f}")

    if args.plot:
        if plt is None:
            print(
                "error: plotting requested but matplotlib is not installed. "
                "Install it with: pip install matplotlib",
                file=sys.stderr,
            )
            return 2

        means = [r[2] for r in rows]
        variances = [r[3] for r in rows]
        sizes = [max(20.0, min(250.0, 15.0 + r[1] * 0.05)) for r in rows]

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(means, variances, s=sizes, alpha=0.75)
        for cls_id, _, mean, var in rows:
            ax.annotate(str(cls_id), (mean, var), fontsize=8, xytext=(4, 3), textcoords="offset points")

        ax.set_title("2D Class Defocus Statistics")
        ax.set_xlabel("Mean Defocus")
        ax.set_ylabel("Variance of Defocus")
        ax.grid(True, linestyle="--", alpha=0.4)
        fig.tight_layout()
        fig.savefig(args.plot_out, dpi=150)
        print(f"Saved plot: {args.plot_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
