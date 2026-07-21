#!/usr/bin/env python3
"""Measure SIMPLE abinitio2D memory use with a controlled screening design.

Each case runs in a fresh ``simple_exec`` process. Preparation creates a
deterministic synthetic particle stack and imports it into a minimal SIMPLE
project before memory telemetry starts.

The default screening design combines one-factor sweeps with a small set of
jointly varied cases. This captures useful interactions without constructing a
potentially enormous full Cartesian product.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import re
import shutil
import struct
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from benchmark_motion_correct import MIB, locate_executable, read_telemetry, run_logged


RESULT_FIELDS = (
    "commander",
    "case_id",
    "design_group",
    "repeat",
    "nptcls",
    "box",
    "particle_pixels",
    "total_input_pixels",
    "smpd_angstrom_per_pixel",
    "physical_box_angstrom",
    "effective_box",
    "effective_particle_pixels",
    "autoscale",
    "nthr",
    "ncls",
    "nrefs",
    "particles_per_reference",
    "nsample_requested",
    "nsample_effective",
    "nstages",
    "nits_per_stage",
    "extr_lim",
    "eo_stage",
    "refine",
    "sigma_est",
    "ctf",
    "mskdiam_angstrom",
    "synthetic_noise_sigma",
    "start_current_rss_bytes",
    "peak_rss_bytes",
    "peak_delta_rss_bytes",
    "peak_rss_mib",
    "peak_delta_rss_mib",
    "delta_bytes_per_input_pixel",
    "elapsed_s",
    "status",
    "returncode",
    "telemetry_file",
    "log_file",
)


@dataclass(frozen=True)
class Case:
    nptcls: int
    box: int
    nthr: int
    ncls: int
    nsample: int = 0
    autoscale: bool = False
    group: str = "screening"

    @property
    def key(self) -> tuple[int, int, int, int, int, bool]:
        return self.nptcls, self.box, self.nthr, self.ncls, self.nsample, self.autoscale

    def case_id(self, repeat: int) -> str:
        sample = "all" if self.nsample == 0 else f"{self.nsample:06d}"
        autoscale = "y" if self.autoscale else "n"
        return (
            f"n{self.nptcls:06d}_b{self.box:04d}_t{self.nthr:02d}_"
            f"c{self.ncls:03d}_s{sample}_a{autoscale}_r{repeat:02d}"
        )


def positive_int(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return value


def nonnegative_int(text: str) -> int:
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError("value must be nonnegative")
    return value


def positive_float(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("value must be finite and positive")
    return value


def unique_cases(cases: Iterable[Case]) -> list[Case]:
    seen: set[tuple[int, int, int, int, int, bool]] = set()
    result: list[Case] = []
    for case in cases:
        if case.key not in seen:
            seen.add(case.key)
            result.append(case)
    return result


def screening_design(args: argparse.Namespace) -> list[Case]:
    baseline = Case(
        args.baseline_particles,
        args.baseline_box,
        args.baseline_threads,
        args.baseline_references,
        group="baseline",
    )
    cases: list[Case] = [baseline]
    cases.extend(
        Case(value, baseline.box, baseline.nthr, baseline.ncls, group="particle_count")
        for value in args.particle_counts
    )
    cases.extend(
        Case(baseline.nptcls, value, baseline.nthr, baseline.ncls, group="box_size")
        for value in args.box_sizes
    )
    cases.extend(
        Case(baseline.nptcls, baseline.box, value, baseline.ncls, group="thread_count")
        for value in args.threads
    )
    cases.extend(
        Case(baseline.nptcls, baseline.box, baseline.nthr, value, group="reference_count")
        for value in args.references
    )

    if args.include_nsample:
        largest_particle_count = max(args.particle_counts)
        cases.extend(
            Case(
                largest_particle_count,
                baseline.box,
                baseline.nthr,
                baseline.ncls,
                nsample=value,
                group="sampled_particles",
            )
            for value in args.nsamples
        )
    cases.extend(
        Case(baseline.nptcls, value, baseline.nthr, baseline.ncls, autoscale=True, group="autoscale")
        for value in sorted({baseline.box, max(args.box_sizes), args.box_sizes[len(args.box_sizes) // 2]})
    )

    if not args.no_interactions:
        interactions = (
            (100, 64, 1, 2, 0),
            (250, 96, 2, 4, 0),
            (500, 128, 4, 8, 0),
            (1000, 160, 8, 16, 0),
            (2000, 192, 4, 32, 0),
            (2000, 64, 8, 16, 0),
            (100, 192, 8, 32, 0),
            (1000, 128, 1, 32, 0),
            (250, 160, 4, 2, 0),
            (500, 192, 2, 16, 0),
            (2000, 128, 2, 8, 0),
            (1000, 96, 8, 4, 0),
        )
        cases.extend(Case(*values, group="interaction") for values in interactions)
    return unique_cases(cases)


def factorial_design(args: argparse.Namespace) -> list[Case]:
    return [
        Case(nptcls, box, nthr, ncls, group="factorial")
        for nptcls in args.particle_counts
        for box in args.box_sizes
        for nthr in args.threads
        for ncls in args.references
    ]


def make_mrc_header(box: int, nptcls: int, smpd: float) -> bytes:
    header = bytearray(1024)
    struct.pack_into("<4i", header, 0, box, box, nptcls, 2)
    struct.pack_into("<3i", header, 28, box, box, nptcls)
    struct.pack_into("<3f", header, 40, box * smpd, box * smpd, nptcls * smpd)
    struct.pack_into("<3f", header, 52, 90.0, 90.0, 90.0)
    struct.pack_into("<3i", header, 64, 1, 2, 3)
    struct.pack_into("<3f", header, 76, -1.0, 1.0, 0.0)
    header[208:212] = b"MAP "
    header[212:216] = bytes((0x44, 0x41, 0x00, 0x00))
    struct.pack_into("<f", header, 216, 0.5)
    struct.pack_into("<i", header, 220, 1)
    label = b"SIMPLE deterministic abinitio2D memory-test particles"
    header[224 : 224 + len(label)] = label
    return bytes(header)


def write_synthetic_particles(
    path: Path,
    nptcls: int,
    box: int,
    smpd: float,
    noise_sigma: float,
) -> None:
    """Write deterministic asymmetric particles in bounded NumPy batches."""
    path.parent.mkdir(parents=True, exist_ok=True)
    axis = np.linspace(-1.0, 1.0, box, dtype=np.float32)
    y, x = np.meshgrid(axis, axis, indexing="ij")
    rng = np.random.default_rng(20260720 + nptcls * 31 + box * 101)
    batch_size = min(64, nptcls)
    with path.open("wb") as stream:
        stream.write(make_mrc_header(box, nptcls, smpd))
        for first in range(0, nptcls, batch_size):
            count = min(batch_size, nptcls - first)
            indices = np.arange(first, first + count, dtype=np.float32)
            angle = indices * np.float32(2.399963229728653)
            cosine = np.cos(angle)[:, None, None]
            sine = np.sin(angle)[:, None, None]
            xr = x[None, :, :] * cosine + y[None, :, :] * sine
            yr = -x[None, :, :] * sine + y[None, :, :] * cosine
            main = np.exp(-(xr * xr / 0.18 + yr * yr / 0.045))
            side = 0.65 * np.exp(-((xr - 0.38) ** 2 / 0.035 + (yr + 0.18) ** 2 / 0.06))
            notch = 0.35 * np.exp(-((xr + 0.30) ** 2 / 0.025 + (yr - 0.25) ** 2 / 0.04))
            texture = noise_sigma * rng.standard_normal((count, box, box), dtype=np.float32)
            particles = main + side - notch + texture
            particles -= particles.mean(axis=(1, 2), keepdims=True)
            scale = particles.std(axis=(1, 2), keepdims=True)
            particles /= np.maximum(scale, np.float32(1.0e-6))
            particles.astype("<f4", copy=False).tofile(stream)


def prepare_case(
    case_dir: Path,
    stack: Path,
    simple_exec: Path,
    smpd: float,
    env: dict[str, str],
) -> Path:
    case_dir.mkdir(parents=True, exist_ok=False)
    prep_log = case_dir / "prepare.log"
    project_name = "abinitio2d_memory"
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
    rc = run_logged(
        [
            str(simple_exec),
            "prg=import_particles",
            f"projfile={project_file}",
            f"stk={stack.resolve()}",
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
        raise RuntimeError(f"import_particles failed with exit status {rc}; see {prep_log}")
    return project_file


def parse_effective_box(log_path: Path, fallback: int) -> int:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    matches = re.findall(r"ORIGINAL/CROPPED IMAGE SIZE \(pixels\):\s*\d+/\s*(\d+)", text)
    return int(matches[-1]) if matches else fallback


def measure_case(
    case_dir: Path,
    project_file: Path,
    simple_exec: Path,
    case: Case,
    repeat: int,
    args: argparse.Namespace,
    output_dir: Path,
    env: dict[str, str],
) -> dict[str, object]:
    measure_dir = case_dir / "measure"
    measure_dir.mkdir()
    log_path = measure_dir / "abinitio2D.log"
    command = [
        str(simple_exec),
        "prg=abinitio2D",
        f"projfile={project_file}",
        "mkdir=no",
        f"ncls={case.ncls}",
        f"mskdiam={args.mskdiam:.8g}",
        f"smpd={args.smpd:.8g}",
        f"nthr={case.nthr}",
        f"nstages={args.nstages}",
        f"nits_per_stage={args.nits_per_stage}",
        f"extr_lim={args.extr_lim}",
        f"eo_stage={args.eo_stage}",
        f"autoscale={'yes' if case.autoscale else 'no'}",
        f"refine={args.refine}",
        f"sigma_est={args.sigma_est}",
        "rank_cavgs=no",
        "memreport=yes",
        "memreport_interval=1",
    ]
    if case.nsample > 0:
        command.append(f"nsample={case.nsample}")
    started = time.monotonic()
    rc = run_logged(command, measure_dir, log_path, env)
    wall_elapsed = time.monotonic() - started
    telemetry_rel = ""
    start_current = peak = delta = 0
    measured_elapsed = wall_elapsed
    try:
        telemetry, start_current, peak, measured_elapsed = read_telemetry(measure_dir)
        telemetry_rel = str(telemetry.relative_to(output_dir))
        delta = max(0, peak - start_current)
    except RuntimeError:
        if rc == 0:
            raise
    effective_box = parse_effective_box(log_path, case.box)
    effective_sample = min(case.nptcls, case.nsample if case.nsample > 0 else case.nptcls)
    total_input_pixels = case.nptcls * case.box * case.box
    return {
        "commander": "abinitio2D",
        "case_id": case.case_id(repeat),
        "design_group": case.group,
        "repeat": repeat,
        "nptcls": case.nptcls,
        "box": case.box,
        "particle_pixels": case.box * case.box,
        "total_input_pixels": total_input_pixels,
        "smpd_angstrom_per_pixel": f"{args.smpd:.8g}",
        "physical_box_angstrom": f"{case.box * args.smpd:.6f}",
        "effective_box": effective_box,
        "effective_particle_pixels": effective_box * effective_box,
        "autoscale": "yes" if case.autoscale else "no",
        "nthr": case.nthr,
        "ncls": case.ncls,
        "nrefs": case.ncls,
        "particles_per_reference": f"{case.nptcls / case.ncls:.6f}",
        "nsample_requested": case.nsample,
        "nsample_effective": effective_sample,
        "nstages": args.nstages,
        "nits_per_stage": args.nits_per_stage,
        "extr_lim": args.extr_lim,
        "eo_stage": args.eo_stage,
        "refine": args.refine,
        "sigma_est": args.sigma_est,
        "ctf": "no",
        "mskdiam_angstrom": f"{args.mskdiam:.6f}",
        "synthetic_noise_sigma": f"{args.noise_sigma:.6f}",
        "start_current_rss_bytes": start_current,
        "peak_rss_bytes": peak,
        "peak_delta_rss_bytes": delta,
        "peak_rss_mib": f"{peak / MIB:.6f}",
        "peak_delta_rss_mib": f"{delta / MIB:.6f}",
        "delta_bytes_per_input_pixel": f"{delta / total_input_pixels:.9f}",
        "elapsed_s": f"{measured_elapsed:.6f}",
        "status": "ok" if rc == 0 else "failed",
        "returncode": rc,
        "telemetry_file": telemetry_rel,
        "log_file": str(log_path.relative_to(output_dir)),
    }


def append_result(path: Path, row: dict[str, object]) -> None:
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RESULT_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def discard_case_data(case_dir: Path) -> None:
    project_dir = case_dir / "abinitio2d_memory"
    if project_dir.exists():
        shutil.rmtree(project_dir)
    measure_dir = case_dir / "measure"
    if measure_dir.exists():
        for path in measure_dir.iterdir():
            if path.name == "abinitio2D.log" or path.name.startswith("memory_usage_"):
                continue
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Collect isolated peak-RSS measurements for SIMPLE abinitio2D.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--design", choices=("screening", "factorial"), default="screening")
    parser.add_argument("--particle-counts", nargs="+", type=positive_int, default=[100, 250, 500, 1000, 2000])
    parser.add_argument("--box-sizes", nargs="+", type=positive_int, default=[64, 96, 128, 160, 192])
    parser.add_argument("--threads", nargs="+", type=positive_int, default=[1, 2, 4, 8])
    parser.add_argument("--references", nargs="+", type=positive_int, default=[2, 4, 8, 16, 32])
    parser.add_argument("--nsamples", nargs="+", type=nonnegative_int, default=[100, 250, 500, 1000, 0])
    parser.add_argument(
        "--include-nsample",
        action="store_true",
        help="include nsample cases; requires at least three stages for a bounded representative workflow",
    )
    parser.add_argument("--baseline-particles", type=positive_int, default=500)
    parser.add_argument("--baseline-box", type=positive_int, default=96)
    parser.add_argument("--baseline-threads", type=positive_int, default=4)
    parser.add_argument("--baseline-references", type=positive_int, default=8)
    parser.add_argument("--smpd", type=positive_float, default=1.3)
    parser.add_argument("--mskdiam", type=positive_float, default=60.0)
    parser.add_argument("--noise-sigma", type=positive_float, default=0.15)
    parser.add_argument("--nstages", type=positive_int, default=1)
    parser.add_argument("--nits-per-stage", type=positive_int, default=1)
    parser.add_argument("--extr-lim", type=positive_int, default=4)
    parser.add_argument("--eo-stage", choices=("yes", "no"), default="no")
    parser.add_argument("--refine", choices=("prob_snhc", "prob", "snhc_smpl"), default="prob_snhc")
    parser.add_argument("--sigma-est", choices=("global", "group"), default="global")
    parser.add_argument("--repeats", type=positive_int, default=1)
    parser.add_argument("--max-cases", type=positive_int, default=100)
    parser.add_argument("--no-interactions", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--discard-case-data", action="store_true")
    parser.add_argument("--discard-inputs", action="store_true")
    parser.add_argument("--simple-exec", help="path to simple_exec")
    parser.add_argument("--output-dir", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    cases = screening_design(args) if args.design == "screening" else factorial_design(args)
    if args.include_nsample and args.nstages < 3:
        parser.error("--include-nsample requires --nstages of at least 3")
    if len(cases) * args.repeats > args.max_cases:
        parser.error(
            f"design contains {len(cases) * args.repeats} runs, exceeding --max-cases={args.max_cases}"
        )
    for case in cases:
        if case.ncls > case.nptcls:
            parser.error(f"ncls={case.ncls} exceeds nptcls={case.nptcls} in {case}")
        if args.mskdiam >= case.box * args.smpd:
            parser.error(
                f"mskdiam={args.mskdiam:g} A must be smaller than the physical box "
                f"({case.box * args.smpd:g} A) for box={case.box}"
            )
    if args.dry_run:
        writer = csv.writer(sys.stdout)
        writer.writerow(("case", "group", "nptcls", "box", "nthr", "ncls", "nsample", "autoscale"))
        for index, case in enumerate(cases, 1):
            writer.writerow((index, case.group, case.nptcls, case.box, case.nthr, case.ncls, case.nsample, case.autoscale))
        return 0

    script_path = Path(__file__)
    simple_exec = locate_executable(args.simple_exec, "simple_exec", script_path)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (args.output_dir or Path(f"abinitio2d_memory_{timestamp}")).expanduser().resolve()
    if output_dir.exists() and not args.resume:
        raise FileExistsError(f"refusing to overwrite existing output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=args.resume)
    cases_dir = output_dir / "cases"
    inputs_dir = output_dir / "inputs"
    cases_dir.mkdir(exist_ok=True)
    inputs_dir.mkdir(exist_ok=True)
    results_path = output_dir / "results.csv"

    completed: set[str] = set()
    if args.resume and results_path.exists():
        with results_path.open(newline="", encoding="utf-8") as stream:
            completed = {row["case_id"] for row in csv.DictReader(stream)}
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "commander": "abinitio2D",
        "design": args.design,
        "planned_unique_cases": len(cases),
        "repeats": args.repeats,
        "simple_exec": str(simple_exec),
        "platform": platform.platform(),
        "python": sys.version,
        "measurement": "peak_rss_bytes minus current_rss_bytes at telemetry start",
        "fixed_configuration": {
            "smpd": args.smpd,
            "mskdiam": args.mskdiam,
            "nstages": args.nstages,
            "nits_per_stage": args.nits_per_stage,
            "extr_lim": args.extr_lim,
            "eo_stage": args.eo_stage,
            "refine": args.refine,
            "sigma_est": args.sigma_est,
            "ctf": "no",
        },
    }
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    env = os.environ.copy()
    failures = 0
    total = len(cases) * args.repeats
    index = 0
    for case in cases:
        stack = inputs_dir / f"particles_n{case.nptcls:06d}_b{case.box:04d}.mrc"
        if not stack.exists():
            write_synthetic_particles(stack, case.nptcls, case.box, args.smpd, args.noise_sigma)
        for repeat in range(1, args.repeats + 1):
            index += 1
            case_id = case.case_id(repeat)
            if case_id in completed:
                print(f"[{index}/{total}] {case_id} skipped (already recorded)", flush=True)
                continue
            case_dir = cases_dir / case_id
            print(f"[{index}/{total}] {case_id} ({case.group})", flush=True)
            env["OMP_NUM_THREADS"] = str(case.nthr)
            try:
                project_file = prepare_case(case_dir, stack, simple_exec, args.smpd, env)
                row = measure_case(case_dir, project_file, simple_exec, case, repeat, args, output_dir, env)
            except Exception as exc:
                failures += 1
                row = {field: "" for field in RESULT_FIELDS}
                row.update(
                    commander="abinitio2D",
                    case_id=case_id,
                    design_group=case.group,
                    repeat=repeat,
                    nptcls=case.nptcls,
                    box=case.box,
                    nthr=case.nthr,
                    ncls=case.ncls,
                    nrefs=case.ncls,
                    nsample_requested=case.nsample,
                    autoscale="yes" if case.autoscale else "no",
                    status="setup_failed",
                    returncode=-1,
                    log_file=str((case_dir / "prepare.log").relative_to(output_dir)),
                )
                case_dir.mkdir(parents=True, exist_ok=True)
                (case_dir / "harness_error.txt").write_text(str(exc) + "\n", encoding="utf-8")
                print(f"  failed: {exc}", file=sys.stderr, flush=True)
            else:
                if row["status"] != "ok":
                    failures += 1
                print(
                    f"  peak={row['peak_rss_mib']} MiB delta={row['peak_delta_rss_mib']} MiB "
                    f"effective_box={row['effective_box']} status={row['status']}",
                    flush=True,
                )
            append_result(results_path, row)
            if args.discard_case_data:
                discard_case_data(case_dir)

    if args.discard_inputs and inputs_dir.exists():
        shutil.rmtree(inputs_dir)
    print(f"Results: {results_path}")
    if failures:
        print(f"{failures} of {total} cases failed; inspect per-case logs", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
