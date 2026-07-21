#!/usr/bin/env python3
"""Measure SIMPLE commander memory use over an input-parameter grid.

The first supported commander is ``motion_correct``.  Each grid point is run
through ``simple_private_exec`` in a fresh process because peak RSS is a
process-lifetime statistic and cannot be reset reliably between in-process
trials.

The harness creates a deterministic, correlated synthetic MRC movie for each
case, imports it into a minimal SIMPLE project, enables SIMPLE's native memory
telemetry, and appends one row to ``results.csv``.  Preparation is deliberately
outside the measured process.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import shutil
import struct
import subprocess
import sys
import time
from array import array
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence


MIB = 1024 * 1024
RESULT_FIELDS = (
    "commander",
    "case_id",
    "repeat",
    "xdim",
    "ydim",
    "pixels_per_frame",
    "effective_xdim",
    "effective_ydim",
    "effective_pixels_per_frame",
    "frames",
    "smpd_angstrom_per_pixel",
    "smpd_downscale_angstrom_per_pixel",
    "physical_width_angstrom",
    "physical_height_angstrom",
    "nthr",
    "start_current_rss_bytes",
    "peak_rss_bytes",
    "peak_delta_rss_bytes",
    "peak_rss_mib",
    "peak_delta_rss_mib",
    "delta_bytes_per_frame_pixel",
    "elapsed_s",
    "status",
    "returncode",
    "telemetry_file",
    "log_file",
)


@dataclass(frozen=True)
class Dimensions:
    x: int
    y: int

    @property
    def label(self) -> str:
        return f"{self.x}x{self.y}"


def positive_int(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return value


def positive_float(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or value <= 0.0:
        raise argparse.ArgumentTypeError("value must be finite and greater than 0")
    return value


def parse_dimensions(text: str) -> Dimensions:
    normalized = text.lower().replace("*", "x")
    if "x" in normalized:
        parts = normalized.split("x")
        if len(parts) != 2:
            raise argparse.ArgumentTypeError(f"invalid dimensions: {text}")
        xdim, ydim = (positive_int(part) for part in parts)
    else:
        xdim = ydim = positive_int(normalized)
    return Dimensions(xdim, ydim)


def locate_executable(explicit: str | None, name: str, script_path: Path) -> Path:
    if explicit:
        candidate = Path(explicit).expanduser().resolve()
    else:
        from_path = shutil.which(name)
        if from_path:
            candidate = Path(from_path).resolve()
        else:
            repo_root = script_path.resolve().parent.parent.parent
            candidates = (
                repo_root / "build" / "bin" / name,
                repo_root / "build" / "production" / name,
            )
            candidate = next((path for path in candidates if path.is_file()), candidates[0])
    if not candidate.is_file() or not os.access(candidate, os.X_OK):
        raise FileNotFoundError(f"executable not found or not executable: {candidate}")
    return candidate


def make_mrc_header(xdim: int, ydim: int, nframes: int, smpd: float) -> bytes:
    """Return a little-endian MRC2014 header for a float32 image stack."""
    header = bytearray(1024)
    struct.pack_into("<4i", header, 0, xdim, ydim, nframes, 2)
    struct.pack_into("<3i", header, 16, 0, 0, 0)
    struct.pack_into("<3i", header, 28, xdim, ydim, nframes)
    struct.pack_into("<3f", header, 40, xdim * smpd, ydim * smpd, nframes * smpd)
    struct.pack_into("<3f", header, 52, 90.0, 90.0, 90.0)
    struct.pack_into("<3i", header, 64, 1, 2, 3)
    struct.pack_into("<3f", header, 76, -1.0, 1.0, 0.0)
    struct.pack_into("<2i", header, 88, 0, 0)
    struct.pack_into("<3f", header, 196, 0.0, 0.0, 0.0)
    header[208:212] = b"MAP "
    header[212:216] = bytes((0x44, 0x41, 0x00, 0x00))
    struct.pack_into("<f", header, 216, 0.45)
    struct.pack_into("<i", header, 220, 1)
    label = b"SIMPLE deterministic memory-test movie"
    header[224 : 224 + len(label)] = label
    return bytes(header)


def synthetic_row(xdim: int, ydim: int, y: int) -> bytes:
    """Create a feature-rich row; all frames share the same image exactly."""
    values = array("f")
    yphase = 2.0 * math.pi * y / max(ydim, 1)
    for x in range(xdim):
        xphase = 2.0 * math.pi * x / max(xdim, 1)
        texture = ((x * 17 + y * 31 + (x ^ y) * 3) % 257) / 256.0 - 0.5
        value = 0.55 * math.sin(7.0 * xphase) + 0.35 * math.cos(5.0 * yphase) + 0.2 * texture
        values.append(value)
    if sys.byteorder != "little":
        values.byteswap()
    return values.tobytes()


def write_synthetic_movie(path: Path, dims: Dimensions, nframes: int, smpd: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame_path = path.with_suffix(path.suffix + ".frame.tmp")
    try:
        with frame_path.open("wb") as frame_stream:
            for y in range(dims.y):
                frame_stream.write(synthetic_row(dims.x, dims.y, y))
        with path.open("wb") as movie_stream:
            movie_stream.write(make_mrc_header(dims.x, dims.y, nframes, smpd))
            for _ in range(nframes):
                with frame_path.open("rb") as frame_stream:
                    shutil.copyfileobj(frame_stream, movie_stream, length=8 * MIB)
    finally:
        if frame_path.exists():
            frame_path.unlink()


def round_to_even(value: float) -> int:
    """Match SIMPLE's round2even for the positive dimensions used here."""
    rounded = math.floor(value + 0.5)
    if rounded % 2 == 0:
        return rounded
    lower = rounded - 1
    upper = rounded + 1
    return lower if abs(lower - value) <= abs(upper - value) else upper


def effective_dimensions(dims: Dimensions, smpd: float, smpd_downscale: float) -> Dimensions:
    scale = min(1.0, smpd / smpd_downscale)
    if scale >= 0.99:
        return dims
    return Dimensions(round_to_even(dims.x * scale), round_to_even(dims.y * scale))


def run_logged(command: Sequence[str], cwd: Path, log_path: Path, env: dict[str, str]) -> int:
    with log_path.open("a", encoding="utf-8") as log:
        log.write("$ " + " ".join(command) + "\n")
        log.flush()
        completed = subprocess.run(
            command,
            cwd=cwd,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
        log.write(f"[exit {completed.returncode}]\n")
    return completed.returncode


def prepare_motion_correct_case(
    case_dir: Path,
    simple_exec: Path,
    dims: Dimensions,
    smpd: float,
    frames: int,
    env: dict[str, str],
) -> Path:
    case_dir.mkdir(parents=True)
    prep_log = case_dir / "prepare.log"
    movie = (case_dir / "synthetic_movie.mrc").resolve()
    write_synthetic_movie(movie, dims, frames, smpd)
    project_name = "motion_correct_memory"
    rc = run_logged(
        [str(simple_exec), "prg=new_project", f"projname={project_name}", "qsys_name=local"],
        case_dir,
        prep_log,
        env,
    )
    if rc != 0:
        raise RuntimeError(f"new_project failed with exit status {rc}; see {prep_log}")
    project_dir = case_dir / project_name
    project_file = (project_dir / f"{project_name}.simple").resolve()
    filetab = project_dir / "movies.txt"
    filetab.write_text(str(movie) + "\n", encoding="utf-8")
    rc = run_logged(
        [
            str(simple_exec),
            "prg=import_movies",
            f"projfile={project_file}",
            f"filetab={filetab.resolve()}",
            f"smpd={smpd:.8g}",
            "kv=300",
            "cs=2.7",
            "fraca=0.1",
            "ctf=no",
            "mkdir=no",
        ],
        project_dir,
        prep_log,
        env,
    )
    if rc != 0:
        raise RuntimeError(f"import_movies failed with exit status {rc}; see {prep_log}")
    return project_file


def read_telemetry(measure_dir: Path) -> tuple[Path, int, int, float]:
    telemetry_files = sorted(measure_dir.glob("memory_usage_*.csv"))
    if len(telemetry_files) != 1:
        raise RuntimeError(f"expected one telemetry CSV in {measure_dir}, found {len(telemetry_files)}")
    telemetry = telemetry_files[0]
    with telemetry.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"telemetry CSV contains no samples: {telemetry}")
    start = next((row for row in rows if row["phase"] == "start"), rows[0])
    start_current = int(start["current_rss_bytes"])
    peak = max(int(row["peak_rss_bytes"]) for row in rows)
    elapsed = max(float(row["elapsed_s"]) for row in rows)
    return telemetry, start_current, peak, elapsed


def append_result(path: Path, row: dict[str, object]) -> None:
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RESULT_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def discard_case_data(case_dir: Path) -> None:
    """Keep measurement evidence while removing large generated case artifacts."""
    movie = case_dir / "synthetic_movie.mrc"
    if movie.exists():
        movie.unlink()
    project_dir = case_dir / "motion_correct_memory"
    if project_dir.exists():
        shutil.rmtree(project_dir)
    measure_dir = case_dir / "measure"
    if measure_dir.exists():
        for path in measure_dir.iterdir():
            if path.name == "motion_correct.log" or path.name.startswith("memory_usage_"):
                continue
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()


def measure_motion_correct(
    case_dir: Path,
    project_file: Path,
    private_exec: Path,
    dims: Dimensions,
    smpd: float,
    frames: int,
    nthr: int,
    smpd_downscale: float,
    repeat: int,
    env: dict[str, str],
) -> dict[str, object]:
    measure_dir = case_dir / "measure"
    measure_dir.mkdir()
    log_path = measure_dir / "motion_correct.log"
    command = [
        str(private_exec),
        "prg=motion_correct",
        f"projfile={project_file}",
        "mkdir=no",
        "algorithm=iso",
        "mcpatch=no",
        "downscale=yes",
        f"smpd_downscale={smpd_downscale:.8g}",
        "fromp=1",
        "top=1",
        "part=1",
        "nparts=1",
        f"nthr={nthr}",
        "outfile=motion_correct_oris.simple",
        "memreport=yes",
        "memreport_interval=1",
    ]
    started = time.monotonic()
    rc = run_logged(command, measure_dir, log_path, env)
    wall_elapsed = time.monotonic() - started
    status = "ok" if rc == 0 else "failed"
    telemetry_rel = ""
    start_current = peak = delta = 0
    measured_elapsed = wall_elapsed
    try:
        telemetry, start_current, peak, measured_elapsed = read_telemetry(measure_dir)
        telemetry_rel = str(telemetry.relative_to(case_dir.parent.parent))
        delta = max(0, peak - start_current)
    except RuntimeError:
        if rc == 0:
            raise
    effective_dims = effective_dimensions(dims, smpd, smpd_downscale)
    frame_pixels = effective_dims.x * effective_dims.y * frames
    return {
        "commander": "motion_correct",
        "case_id": case_dir.name,
        "repeat": repeat,
        "xdim": dims.x,
        "ydim": dims.y,
        "pixels_per_frame": dims.x * dims.y,
        "effective_xdim": effective_dims.x,
        "effective_ydim": effective_dims.y,
        "effective_pixels_per_frame": effective_dims.x * effective_dims.y,
        "frames": frames,
        "smpd_angstrom_per_pixel": f"{smpd:.8g}",
        "smpd_downscale_angstrom_per_pixel": f"{smpd_downscale:.8g}",
        "physical_width_angstrom": f"{dims.x * smpd:.6f}",
        "physical_height_angstrom": f"{dims.y * smpd:.6f}",
        "nthr": nthr,
        "start_current_rss_bytes": start_current,
        "peak_rss_bytes": peak,
        "peak_delta_rss_bytes": delta,
        "peak_rss_mib": f"{peak / MIB:.6f}",
        "peak_delta_rss_mib": f"{delta / MIB:.6f}",
        "delta_bytes_per_frame_pixel": f"{delta / frame_pixels:.9f}",
        "elapsed_s": f"{measured_elapsed:.6f}",
        "status": status,
        "returncode": rc,
        "telemetry_file": telemetry_rel,
        "log_file": str(log_path.relative_to(case_dir.parent.parent)),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Collect isolated peak-RSS measurements for a SIMPLE commander.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--commander", choices=("motion_correct",), default="motion_correct")
    parser.add_argument(
        "--dimensions",
        nargs="+",
        type=parse_dimensions,
        default=[Dimensions(1024, 1024), Dimensions(2048, 2048)],
        metavar="N|NxM",
        help="movie dimensions in pixels",
    )
    parser.add_argument(
        "--smpds",
        nargs="+",
        type=positive_float,
        default=[0.8, 1.2, 2.0],
        metavar="ANGSTROM_PER_PIXEL",
    )
    parser.add_argument(
        "--frames", nargs="+", type=positive_int, default=[8], metavar="COUNT",
        help="movie frame counts (crossed with dimensions, pixel sizes, and threads)",
    )
    parser.add_argument(
        "--threads", nargs="+", type=positive_int, default=[1], metavar="COUNT",
        help="OpenMP thread counts (crossed with dimensions, pixel sizes, and frames)",
    )
    parser.add_argument("--repeats", type=positive_int, default=1)
    parser.add_argument(
        "--smpd-downscale",
        type=positive_float,
        default=1.3,
        help="target sampling distance used by motion_correct movie downscaling",
    )
    parser.add_argument("--simple-exec", help="path to simple_exec")
    parser.add_argument("--private-exec", help="path to simple_private_exec")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="new directory for inputs, logs, telemetry, and results.csv",
    )
    parser.add_argument(
        "--discard-case-data",
        action="store_true",
        help="after each case, retain logs/telemetry but remove generated movies and projects",
    )
    parser.add_argument(
        "--skip-invalid",
        action="store_true",
        help="record unsupported grid points instead of rejecting the complete grid",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    for dims in args.dimensions:
        for smpd in args.smpds:
            effective = effective_dimensions(dims, smpd, args.smpd_downscale)
            if min(effective.x, effective.y) < 512 and not args.skip_invalid:
                parser.error(
                    f"{dims.label} at smpd={smpd:g} becomes {effective.label} after "
                    f"downscaling to {args.smpd_downscale:g} A/pixel; motion_correct "
                    "requires effective dimensions of at least 512x512"
                )
    script_path = Path(__file__)
    simple_exec = locate_executable(args.simple_exec, "simple_exec", script_path)
    private_exec = locate_executable(args.private_exec, "simple_private_exec", script_path)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (args.output_dir or Path(f"motion_correct_memory_{timestamp}")).expanduser().resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {output_dir}")
    output_dir.mkdir(parents=True)
    cases_dir = output_dir / "cases"
    cases_dir.mkdir()
    results_path = output_dir / "results.csv"
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "commander": args.commander,
        "dimensions": [dims.label for dims in args.dimensions],
        "smpds_angstrom_per_pixel": args.smpds,
        "smpd_downscale_angstrom_per_pixel": args.smpd_downscale,
        "frames": args.frames,
        "threads": args.threads,
        "repeats": args.repeats,
        "simple_exec": str(simple_exec),
        "simple_private_exec": str(private_exec),
        "platform": platform.platform(),
        "python": sys.version,
        "measurement": "peak_rss_bytes minus current_rss_bytes at telemetry start",
        "motion_correct_configuration": {
            "algorithm": "iso",
            "mcpatch": "no",
            "downscale": "yes",
        },
    }
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    env = os.environ.copy()
    failures = 0
    total = (len(args.dimensions) * len(args.smpds) * len(args.frames)
             * len(args.threads) * args.repeats)
    index = 0
    for dims in args.dimensions:
        for smpd in args.smpds:
            for frames in args.frames:
                for nthr in args.threads:
                    env["OMP_NUM_THREADS"] = str(nthr)
                    for repeat in range(1, args.repeats + 1):
                        index += 1
                        smpd_label = f"{smpd:.6g}".replace(".", "p")
                        case_id = (f"{dims.label}_smpd{smpd_label}_f{frames:03d}_"
                                   f"t{nthr:02d}_r{repeat:02d}")
                        failures += run_motion_correct_case(
                            args, dims, smpd, frames, nthr, repeat, case_id,
                            index, total, cases_dir, results_path, simple_exec,
                            private_exec, env,
                        )
    print(f"Results: {results_path}")
    if failures:
        print(f"{failures} of {total} cases failed; inspect per-case logs", file=sys.stderr)
        return 1
    return 0


def run_motion_correct_case(
    args: argparse.Namespace,
    dims: Dimensions,
    smpd: float,
    frames: int,
    nthr: int,
    repeat: int,
    case_id: str,
    index: int,
    total: int,
    cases_dir: Path,
    results_path: Path,
    simple_exec: Path,
    private_exec: Path,
    env: dict[str, str],
) -> int:
    """Prepare, measure, and record one point in the parameter grid."""
    case_dir = cases_dir / case_id
    output_dir = cases_dir.parent
    print(f"[{index}/{total}] {case_id}", flush=True)
    effective = effective_dimensions(dims, smpd, args.smpd_downscale)
    if min(effective.x, effective.y) < 512:
        row = {field: "" for field in RESULT_FIELDS}
        row.update(
            commander=args.commander, case_id=case_id, repeat=repeat,
            xdim=dims.x, ydim=dims.y, pixels_per_frame=dims.x * dims.y,
            effective_xdim=effective.x, effective_ydim=effective.y,
            effective_pixels_per_frame=effective.x * effective.y, frames=frames,
            smpd_angstrom_per_pixel=f"{smpd:.8g}",
            smpd_downscale_angstrom_per_pixel=f"{args.smpd_downscale:.8g}",
            physical_width_angstrom=f"{dims.x * smpd:.6f}",
            physical_height_angstrom=f"{dims.y * smpd:.6f}", nthr=nthr,
            status="unsupported_effective_dimension",
        )
        append_result(results_path, row)
        print(f"  skipped: effective dimensions {effective.label} are below 512x512", flush=True)
        return 0
    try:
        project_file = prepare_motion_correct_case(case_dir, simple_exec, dims, smpd, frames, env)
        row = measure_motion_correct(
            case_dir, project_file, private_exec, dims, smpd, frames, nthr,
            args.smpd_downscale, repeat, env,
        )
    except Exception as exc:  # preserve partial grids and identify the failed case
        row = {field: "" for field in RESULT_FIELDS}
        row.update(
            commander=args.commander, case_id=case_id, repeat=repeat,
            xdim=dims.x, ydim=dims.y, pixels_per_frame=dims.x * dims.y,
            effective_xdim=effective.x, effective_ydim=effective.y,
            effective_pixels_per_frame=effective.x * effective.y, frames=frames,
            smpd_angstrom_per_pixel=f"{smpd:.8g}",
            smpd_downscale_angstrom_per_pixel=f"{args.smpd_downscale:.8g}",
            nthr=nthr, status="setup_failed", returncode=-1,
            log_file=str((case_dir / "prepare.log").relative_to(output_dir)),
        )
        case_dir.mkdir(parents=True, exist_ok=True)
        (case_dir / "harness_error.txt").write_text(str(exc) + "\n", encoding="utf-8")
        print(f"  failed: {exc}", file=sys.stderr, flush=True)
        failed = 1
    else:
        failed = int(row["status"] != "ok")
        print(
            f"  peak={row['peak_rss_mib']} MiB delta={row['peak_delta_rss_mib']} MiB "
            f"status={row['status']}", flush=True,
        )
    append_result(results_path, row)
    if args.discard_case_data:
        discard_case_data(case_dir)
    return failed


if __name__ == "__main__":
    raise SystemExit(main())
