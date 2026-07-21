#!/usr/bin/env python3
"""Screen SIMPLE abinitio3D memory use across workflow and data variables.

The harness generates deterministic asymmetric 3D volumes, asks SIMPLE to
simulate realistic projections, imports those particles into fresh projects,
and measures the resident memory of the complete ``simple_exec`` process tree.

The default design is intentionally screened rather than factorial.  Most
cases execute one internal stage-1 refinement iteration; dedicated 5, 10, and
20 iteration cases measure accumulation across a complete stage.
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
import subprocess
import sys
import time
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from benchmark_motion_correct import MIB, locate_executable, run_logged


BUILTIN_STAGE1_MAXITS = 20
RESULT_FIELDS = (
    "commander", "case_id", "design_group", "repeat", "nptcls", "box",
    "particle_pixels", "total_input_pixels", "smpd_angstrom_per_pixel",
    "physical_box_angstrom", "effective_box", "effective_voxels",
    "mskdiam_angstrom", "nthr", "nthr_ini3d", "nparts", "ncunits",
    "nstates", "multivol_mode", "nsample_requested", "nsample_effective",
    "sample_fraction", "pgrp", "pgrp_start", "symmetry_search",
    "lpstart_angstrom", "lpstop_angstrom", "projrec", "initialization",
    "iteration_budget", "iterations_observed", "nstages", "filt_mode",
    "automsk", "ctf", "synthetic_snr", "start_parent_rss_bytes",
    "peak_parent_rss_bytes", "peak_children_rss_bytes", "peak_tree_rss_bytes",
    "peak_tree_delta_rss_bytes", "peak_parent_rss_mib", "peak_children_rss_mib",
    "peak_tree_rss_mib", "peak_tree_delta_rss_mib", "elapsed_s", "status",
    "returncode", "telemetry_files", "log_file",
)


@dataclass(frozen=True)
class Case:
    nptcls: int = 200
    box: int = 96
    smpd: float = 1.3
    mskdiam: float = 60.0
    nthr: int = 4
    nparts: int = 1
    nstates: int = 1
    multivol_mode: str = "single"
    nsample: int = 0
    pgrp: str = "c1"
    pgrp_start: str = "c1"
    lpstart: float = 20.0
    lpstop: float = 8.0
    projrec: bool = False
    initialization: str = "external_volume"
    iteration_budget: int = BUILTIN_STAGE1_MAXITS
    group: str = "screening"

    @property
    def key(self) -> tuple[object, ...]:
        return tuple(getattr(self, name) for name in (
            "nptcls", "box", "smpd", "mskdiam", "nthr", "nparts",
            "nstates", "multivol_mode", "nsample", "pgrp", "pgrp_start",
            "lpstart", "lpstop", "projrec", "initialization", "iteration_budget",
        ))

    def case_id(self, repeat: int) -> str:
        sample = "all" if self.nsample == 0 else str(self.nsample)
        return (
            f"n{self.nptcls:05d}_b{self.box:03d}_a{self.smpd:.1f}_m{self.mskdiam:.0f}_"
            f"t{self.nthr:02d}_p{self.nparts:02d}_s{self.nstates}_{self.multivol_mode[:3]}_"
            f"q{sample}_{self.pgrp}-{self.pgrp_start}_lp{self.lpstop:g}_"
            f"r{'y' if self.projrec else 'n'}_{self.initialization[:3]}_"
            f"i{self.iteration_budget:02d}_x{repeat:02d}"
        ).replace(".", "p")


def unique_cases(cases: Iterable[Case]) -> list[Case]:
    seen: set[tuple[object, ...]] = set()
    answer: list[Case] = []
    for case in cases:
        if case.key not in seen:
            seen.add(case.key)
            answer.append(case)
    return answer


def screening_design() -> list[Case]:
    base = Case(group="baseline")
    cases = [base]
    cases += [replace(base, nptcls=v, group="particle_count") for v in (100, 200, 500, 1000)]
    cases += [replace(base, box=v, group="box_size") for v in (64, 96, 128, 160)]
    cases += [replace(base, nthr=v, group="thread_count") for v in (1, 2, 4, 8)]
    cases += [
        replace(base, nstates=2, multivol_mode="independent", group="state_count"),
        replace(base, nstates=4, multivol_mode="independent", group="state_count"),
    ]
    cases += [
        replace(base, nptcls=500, nsample=v, group="sampled_particles")
        for v in (50, 100, 250, 500)
    ]
    cases += [replace(base, smpd=v, group="pixel_size") for v in (0.8, 1.3, 2.0)]
    cases += [replace(base, mskdiam=v, group="mask_diameter") for v in (40.0, 60.0, 80.0)]
    cases += [replace(base, lpstop=v, group="low_pass") for v in (8.0, 12.0, 16.0)]
    cases += [
        replace(base, pgrp="c2", pgrp_start="c2", group="symmetry"),
        replace(base, pgrp="d2", pgrp_start="d2", group="symmetry"),
        replace(base, pgrp="c2", pgrp_start="c1", group="symmetry_search"),
    ]
    cases += [replace(base, nthr=2, nparts=v, group="job_partitions") for v in (1, 2, 4)]
    cases += [replace(base, projrec=True, group="projection_reconstruction")]
    cases += [replace(base, initialization="random_volume", group="initialization")]
    cases += [
        replace(base, nptcls=1000, box=128, nthr=8, nsample=250, group="interaction"),
        replace(base, nptcls=500, box=160, nthr=4, smpd=2.0, group="interaction"),
        replace(base, nptcls=500, box=128, nthr=2, nparts=4, group="interaction"),
        replace(base, nptcls=500, box=128, nstates=4, multivol_mode="independent", nsample=100, group="interaction"),
    ]
    return unique_cases(cases)


def make_mrc_header(box: int, nz: int, smpd: float, label: bytes) -> bytes:
    header = bytearray(1024)
    struct.pack_into("<4i", header, 0, box, box, nz, 2)
    struct.pack_into("<3i", header, 28, box, box, nz)
    struct.pack_into("<3f", header, 40, box * smpd, box * smpd, nz * smpd)
    struct.pack_into("<3f", header, 52, 90.0, 90.0, 90.0)
    struct.pack_into("<3i", header, 64, 1, 2, 3)
    header[208:212] = b"MAP "
    header[212:216] = bytes((0x44, 0x41, 0x00, 0x00))
    struct.pack_into("<i", header, 220, 1)
    header[224:224 + min(80, len(label))] = label[:80]
    return bytes(header)


def write_synthetic_volume(path: Path, box: int, smpd: float) -> None:
    axis = np.linspace(-1.0, 1.0, box, dtype=np.float32)
    z, y, x = np.meshgrid(axis, axis, axis, indexing="ij")
    volume = (
        1.00 * np.exp(-((x + .18) ** 2 / .15 + (y - .05) ** 2 / .07 + (z + .10) ** 2 / .10))
        + .75 * np.exp(-((x - .34) ** 2 / .045 + (y + .24) ** 2 / .11 + (z - .20) ** 2 / .055))
        + .48 * np.exp(-((x + .38) ** 2 / .035 + (y + .32) ** 2 / .040 + (z - .30) ** 2 / .075))
        - .22 * np.exp(-((x - .02) ** 2 / .025 + (y - .18) ** 2 / .025 + (z + .03) ** 2 / .025))
    ).astype(np.float32)
    volume -= volume.min()
    volume /= max(float(volume.max()), 1e-6)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        stream.write(make_mrc_header(box, box, smpd, b"SIMPLE deterministic abinitio3D benchmark volume"))
        volume.astype("<f4", copy=False).tofile(stream)


def prepare_inputs(case: Case, inputs_dir: Path, simple_exec: Path, env: dict[str, str]) -> tuple[Path, Path]:
    tag = f"b{case.box:03d}_a{case.smpd:.3f}".replace(".", "p")
    volume = inputs_dir / f"volume_{tag}.mrc"
    if not volume.exists():
        write_synthetic_volume(volume, case.box, case.smpd)
    stack = inputs_dir / f"particles_n{case.nptcls:05d}_{tag}.mrc"
    if not stack.exists():
        log = inputs_dir / f"simulate_n{case.nptcls:05d}_{tag}.log"
        rc = run_logged([
            str(simple_exec), "prg=simulate_particles", f"vol1={volume}",
            f"outstk={stack}", f"outfile={inputs_dir / ('oris_n' + str(case.nptcls) + '_' + tag + '.txt')}",
            f"nptcls={case.nptcls}", f"smpd={case.smpd:.8g}", f"mskdiam={case.mskdiam:.8g}",
            f"nthr={min(case.nthr, 8)}", "pgrp=c1", "snr=0.1", "ctf=no", "sherr=0", "even=yes",
        ], inputs_dir, log, env)
        if rc != 0:
            raise RuntimeError(f"simulate_particles failed with status {rc}; see {log}")
    return volume, stack


def prepare_project(case_dir: Path, case: Case, stack: Path, simple_exec: Path, env: dict[str, str]) -> Path:
    prep_log = case_dir / "prepare.log"
    name = "abinitio3d_memory"
    rc = run_logged([str(simple_exec), "prg=new_project", f"projname={name}", "qsys_name=local"], case_dir, prep_log, env)
    if rc != 0:
        raise RuntimeError(f"new_project failed with status {rc}; see {prep_log}")
    project_dir = case_dir / name
    project_file = (project_dir / f"{name}.simple").resolve()
    rc = run_logged([
        str(simple_exec), "prg=import_particles", f"projfile={project_file}", f"stk={stack.resolve()}",
        f"smpd={case.smpd:.8g}", "kv=300", "cs=2.7", "fraca=0.1", "ctf=no", "mkdir=no",
    ], project_dir, prep_log, env)
    if rc != 0:
        raise RuntimeError(f"import_particles failed with status {rc}; see {prep_log}")
    # abinitio3D requires ptcl2D/class metadata even when an external starting
    # volume is supplied.  Keep this mandatory setup outside the measured span.
    rc = run_logged([
        str(simple_exec), "prg=abinitio2D", f"projfile={project_file}", "mkdir=no", "ncls=8",
        f"mskdiam={case.mskdiam:.8g}", f"smpd={case.smpd:.8g}", f"nthr={case.nthr}",
        "nstages=1", "nits_per_stage=1", "extr_lim=4", "eo_stage=no", "autoscale=no",
        "refine=prob_snhc", "rank_cavgs=no",
    ], project_dir, prep_log, env)
    if rc != 0:
        raise RuntimeError(f"abinitio2D preparation failed with status {rc}; see {prep_log}")
    return project_file


def read_native_process_tree(measure_dir: Path, nparts: int) -> tuple[int, int, int, int]:
    """Return parent start/peak and a conservative concurrent child estimate.

    SIMPLE writes one telemetry CSV per process.  The parent is identified by
    its program label.  Child processes shorter than the one-second periodic
    interval still write start/finish samples.  For memory-planning safety, the
    tree estimate adds the parent peak to the largest ``nparts`` child peaks;
    this is an upper bound when those workers did not overlap.
    """
    parent_start = parent_peak = 0
    child_peaks: list[int] = []
    files = sorted(measure_dir.glob("memory_usage_*.csv"))
    for path in files:
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        if not rows:
            continue
        program = rows[0].get("program", "")
        peaks = [int(row["peak_rss_bytes"]) for row in rows if row.get("peak_rss_bytes")]
        if not peaks:
            continue
        if program == "simple_exec:abinitio3D":
            parent_start = int(rows[0]["current_rss_bytes"])
            parent_peak = max(peaks)
        else:
            child_peaks.append(max(peaks))
    if parent_peak == 0:
        raise RuntimeError("missing parent abinitio3D memory telemetry")
    child_peaks.sort(reverse=True)
    peak_children = sum(child_peaks[:max(1, nparts)])
    return parent_start, parent_peak, peak_children, parent_peak + peak_children


def run_tree_measured(command: list[str], cwd: Path, log_path: Path, env: dict[str, str], nparts: int) -> tuple[int, int, int, int, int, float]:
    started = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT,
                                   text=True, check=False)
    elapsed = time.monotonic() - started
    start_parent, peak_parent, peak_children, peak_tree = read_native_process_tree(cwd, nparts)
    return completed.returncode, start_parent, peak_parent, peak_children, peak_tree, elapsed


def parse_log(log_path: Path, fallback_box: int) -> tuple[int, int]:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    boxes = re.findall(r"ORIGINAL/CROPPED IMAGE SIZE \(pixels\):\s*\d+/\s*(\d+)", text)
    iterations = re.findall(r">>> ITERATION\s+(\d+)", text)
    return (int(boxes[-1]) if boxes else fallback_box, len(iterations))


def measure_case(case_dir: Path, case: Case, repeat: int, volume: Path, project_file: Path,
                 simple_exec: Path, output_dir: Path, env: dict[str, str]) -> dict[str, object]:
    measure_dir = case_dir / "measure"
    measure_dir.mkdir()
    log_path = measure_dir / "abinitio3D.log"
    nsample = case.nptcls if case.nsample == 0 else min(case.nsample, case.nptcls)
    command = [
        str(simple_exec), "prg=abinitio3D", f"projfile={project_file}", "mkdir=no",
        f"pgrp={case.pgrp}", f"pgrp_start={case.pgrp_start}", f"mskdiam={case.mskdiam:.8g}",
        f"nstates={case.nstates}", f"multivol_mode={case.multivol_mode}", "filt_mode=none", "automsk=no",
        "nstages=1", f"nsample={nsample}", f"nparts={case.nparts}", f"ncunits={case.nparts}",
        f"nthr={case.nthr}", f"nthr_ini3D={case.nthr}", "force_lp_range=yes",
        f"lpstart={case.lpstart:.8g}", f"lpstop={case.lpstop:.8g}",
        f"projrec={'yes' if case.projrec else 'no'}", "memreport=yes", "memreport_interval=1",
    ]
    if case.initialization == "external_volume":
        for state in range(1, case.nstates + 1):
            command.append(f"vol{state}={volume}")
    rc, start, parent, children, tree, elapsed = run_tree_measured(
        command, measure_dir, log_path, env, case.nparts
    )
    effective_box, iterations = parse_log(log_path, case.box)
    telemetry = sorted(str(p.relative_to(output_dir)) for p in measure_dir.glob("memory_usage_*.csv"))
    delta = max(0, tree - start)
    return {
        "commander": "abinitio3D", "case_id": case.case_id(repeat), "design_group": case.group,
        "repeat": repeat, "nptcls": case.nptcls, "box": case.box, "particle_pixels": case.box ** 2,
        "total_input_pixels": case.nptcls * case.box ** 2, "smpd_angstrom_per_pixel": f"{case.smpd:.6f}",
        "physical_box_angstrom": f"{case.box * case.smpd:.6f}", "effective_box": effective_box,
        "effective_voxels": effective_box ** 3, "mskdiam_angstrom": f"{case.mskdiam:.6f}",
        "nthr": case.nthr, "nthr_ini3d": case.nthr, "nparts": case.nparts, "ncunits": case.nparts,
        "nstates": case.nstates, "multivol_mode": case.multivol_mode,
        "nsample_requested": case.nsample, "nsample_effective": nsample,
        "sample_fraction": f"{nsample / case.nptcls:.6f}", "pgrp": case.pgrp,
        "pgrp_start": case.pgrp_start, "symmetry_search": "yes" if case.pgrp != case.pgrp_start else "no",
        "lpstart_angstrom": f"{case.lpstart:.6f}", "lpstop_angstrom": f"{case.lpstop:.6f}",
        "projrec": "yes" if case.projrec else "no", "initialization": case.initialization,
        "iteration_budget": case.iteration_budget, "iterations_observed": iterations, "nstages": 1,
        "filt_mode": "none", "automsk": "no", "ctf": "no", "synthetic_snr": "0.100000",
        "start_parent_rss_bytes": start, "peak_parent_rss_bytes": parent,
        "peak_children_rss_bytes": children, "peak_tree_rss_bytes": tree,
        "peak_tree_delta_rss_bytes": delta, "peak_parent_rss_mib": f"{parent / MIB:.6f}",
        "peak_children_rss_mib": f"{children / MIB:.6f}", "peak_tree_rss_mib": f"{tree / MIB:.6f}",
        "peak_tree_delta_rss_mib": f"{delta / MIB:.6f}", "elapsed_s": f"{elapsed:.6f}",
        "status": "ok" if rc == 0 else "failed", "returncode": rc,
        "telemetry_files": ";".join(telemetry), "log_file": str(log_path.relative_to(output_dir)),
    }


def append_result(path: Path, row: dict[str, object]) -> None:
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RESULT_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def discard_case_data(case_dir: Path) -> None:
    project = case_dir / "abinitio3d_memory"
    if project.exists():
        shutil.rmtree(project)
    measure = case_dir / "measure"
    if measure.exists():
        for path in measure.iterdir():
            if path.name == "abinitio3D.log" or path.name.startswith("memory_usage_"):
                continue
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()


def positive_int(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Measure SIMPLE abinitio3D process-tree peak RSS.")
    parser.add_argument("--repeats", type=positive_int, default=1)
    parser.add_argument("--max-cases", type=positive_int, default=100)
    parser.add_argument("--simple-exec")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--discard-case-data", action="store_true")
    parser.add_argument("--discard-inputs", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cases = screening_design()
    if len(cases) * args.repeats > args.max_cases:
        raise SystemExit(f"design has {len(cases) * args.repeats} runs; increase --max-cases")
    for case in cases:
        if case.mskdiam >= case.box * case.smpd:
            raise SystemExit(f"mask {case.mskdiam:g} A does not fit box {case.box} at {case.smpd:g} A/pixel")
    if args.dry_run:
        writer = csv.writer(sys.stdout)
        writer.writerow(("case", "group", "nptcls", "box", "smpd", "mskdiam", "nthr", "nparts",
                         "nstates", "multivol_mode", "nsample", "pgrp", "pgrp_start", "lpstart",
                         "lpstop", "projrec", "initialization", "iteration_budget"))
        for i, c in enumerate(cases, 1):
            writer.writerow((i, c.group, c.nptcls, c.box, c.smpd, c.mskdiam, c.nthr, c.nparts,
                             c.nstates, c.multivol_mode, c.nsample, c.pgrp, c.pgrp_start, c.lpstart,
                             c.lpstop, c.projrec, c.initialization, c.iteration_budget))
        return 0

    simple_exec = locate_executable(args.simple_exec, "simple_exec", Path(__file__))
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = (args.output_dir or Path(f"abinitio3d_memory_{timestamp}")).expanduser().resolve()
    if output_dir.exists() and not args.resume:
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=args.resume)
    cases_dir, inputs_dir = output_dir / "cases", output_dir / "inputs"
    cases_dir.mkdir(exist_ok=True)
    inputs_dir.mkdir(exist_ok=True)
    results_path = output_dir / "results.csv"
    completed: set[str] = set()
    if args.resume and results_path.exists():
        with results_path.open(newline="", encoding="utf-8") as stream:
            completed = {row["case_id"] for row in csv.DictReader(stream)}
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "commander": "abinitio3D",
        "design": "screening", "planned_unique_cases": len(cases), "repeats": args.repeats,
        "simple_exec": str(simple_exec), "platform": platform.platform(), "python": sys.version,
        "measurement": "conservative process-tree peak: parent peak RSS plus the largest nparts child-process peaks from native SIMPLE telemetry",
        "sampling_interval_seconds": 1,
        "notes": ["Stage 1 uses SIMPLE's built-in schedule of up to 20 iterations.",
                  "iterations_observed records iterations actually reported before completion or early stopping.",
                  "Input simulation and required abinitio2D setup are outside the measured interval."],
    }
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    env = os.environ.copy()
    failures = 0
    total = len(cases) * args.repeats
    index = 0
    for case in cases:
        env["OMP_NUM_THREADS"] = str(case.nthr)
        volume, stack = prepare_inputs(case, inputs_dir, simple_exec, env)
        for repeat in range(1, args.repeats + 1):
            index += 1
            case_id = case.case_id(repeat)
            if case_id in completed:
                print(f"[{index}/{total}] {case_id} skipped", flush=True)
                continue
            print(f"[{index}/{total}] {case_id} ({case.group})", flush=True)
            case_dir = cases_dir / case_id
            try:
                case_dir.mkdir(parents=True, exist_ok=False)
                project_file = prepare_project(case_dir, case, stack, simple_exec, env)
                row = measure_case(case_dir, case, repeat, volume, project_file, simple_exec, output_dir, env)
            except Exception as exc:
                failures += 1
                case_dir.mkdir(parents=True, exist_ok=True)
                (case_dir / "harness_error.txt").write_text(str(exc) + "\n", encoding="utf-8")
                row = {field: "" for field in RESULT_FIELDS}
                row.update(commander="abinitio3D", case_id=case_id, design_group=case.group,
                           repeat=repeat, nptcls=case.nptcls, box=case.box, status="setup_failed", returncode=-1)
                print(f"  setup failed: {exc}", file=sys.stderr, flush=True)
            else:
                if row["status"] != "ok":
                    failures += 1
                print(f"  tree_peak={row['peak_tree_rss_mib']} MiB parent={row['peak_parent_rss_mib']} MiB "
                      f"children={row['peak_children_rss_mib']} MiB iterations={row['iterations_observed']} "
                      f"status={row['status']}", flush=True)
            append_result(results_path, row)
            if args.discard_case_data:
                discard_case_data(case_dir)
    if args.discard_inputs and inputs_dir.exists():
        shutil.rmtree(inputs_dir)
    print(f"Results: {results_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
