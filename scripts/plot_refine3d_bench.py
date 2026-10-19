#!/usr/bin/env python3
"""Summarize the benchmark reports of a SIMPLE run: tables, a text report and a figure.

Usage:
    python3 scripts/plot_refine3d_bench.py RUN_OR_BENCH_DIRECTORY [OUTPUT_DIRECTORY]

The input is a run directory holding BENCH/ (one file per benchmark family, and
per partition for the partition families, one record appended per report), the
BENCH/ directory itself, or a directory of the per-iteration report files of
earlier builds (*_BENCH_ITERnnn[_PARTppp].txt and SOLVE3D_BENCH.txt). The
output goes to OUTPUT_DIRECTORY, by default bench_summary/ beside the reports:

    bench_records.csv     every value of every record, one row each (long table)
    bench_iterations.csv  one row per refine3D iteration: the strategy phases,
                          their critical-path split, the slowest and mean
                          partition times, the mean percentage of the
                          search space scanned, and the nonuniform filter's
                          split in volume assembly
    bench_report.txt      the run, the time split by stage, the overheads, the
                          load balance, the outlier iterations, the share of the 2D
                          tables in the refine2D iterations, the thread-seconds
                          of the probability tables and the searches, the work
                          of the continuous pose solves (refine=cont, the
                          polish), the pre-step driver per number of states and
                          directions, and the worker phases
    bench_profile.svg     the per-iteration critical path (stacked) and the
                          slowest and mean partition times
    bench_iterations2D.csv, bench_profile2D.svg
                          the same for the refine2D iterations (solve2D's
                          stages), from the REFINE2D_STRATEGY reports

Only the Python standard library is used.

The critical path of an iteration splits the strategy's wall time into what
the slowest partition did and what the master spent around it:

    setup/init                    strategy setup, reprojection model included
    probability tables            the slowest probability-table part
    pre-step outside the tables   pre-step minus the slowest part: launch,
                                  waiting, table reads, ref_assign, writes
    search                        the slowest search worker
    scheduling outside the search matcher/scheduler minus the slowest worker:
                                  launch, the scheduler's polling and the
                                  sigma2 candidate preparation
    volume assembly               the volassemble process
    assembly outside volassemble  assembly/postprocess minus volassemble
    sigma2 commit                 the commit after the assembly
    unaccounted                   the strategy total minus all of the above

The 2D critical path splits a refine2D iteration (to the end of its
finalisation) the same way: setup, the slowest 2D table part, the pre-step
outside it, the slowest search partition, the scheduling outside it (in shared
memory the matcher runs in the master, its class sums and restoration
included), class-average assembly (cavgassemble) and the master's time around
it, the merge of the partitions' documents, the master's project writes and
rereads, the sigma2 preparation and commit, convergence with the star file and
the JPEG, and the unaccounted rest.
"""

from __future__ import annotations

import csv
import math
import re
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from html import escape
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

NUMBER_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$")
HEADER_RE = re.compile(r"^\s*\*{3}\s*(.*?)\s*\*{3}\s*$")
VALUE_RE = re.compile(r"^\s*(.+?)\s*:\s*(.*?)\s*$")
OLD_NAME_RE = re.compile(r"^(?P<family>[A-Z0-9_]+?)_BENCH(?:_ITER(?P<iter>\d+))?(?:_PART(?P<part>\d+))?\.txt$")
NEW_NAME_RE = re.compile(r"^(?P<family>[A-Z0-9_]+?)_BENCH(?:_PART(?P<part>\d+))?\.txt$")


@dataclass
class Record:
    """One benchmark report: its family, iteration, partition and pass, and its sections."""
    family: str
    source: str
    iteration: Optional[int] = None
    partition: Optional[int] = None
    npass: int = 1
    context: Dict[str, str] = field(default_factory=dict)
    seconds: Dict[str, float] = field(default_factory=dict)
    percent: Dict[str, float] = field(default_factory=dict)
    thread_seconds: Dict[str, float] = field(default_factory=dict)
    cont_work: Dict[str, float] = field(default_factory=dict)

    def empty(self) -> bool:
        return not (self.context or self.seconds or self.percent or self.thread_seconds or self.cont_work)


def as_number(text: str) -> Optional[float]:
    text = text.strip()
    return float(text) if NUMBER_RE.match(text) else None


def section_kind(title: str) -> str:
    title = " ".join(title.upper().split())
    if title == "BENCHMARK RECORD":
        return "record"
    if title == "BENCHMARK CONTEXT":
        return "context"
    if title.endswith("TIMINGS (S)"):
        return "seconds"
    if title.startswith("RELATIVE TIMINGS"):
        return "percent"
    if title.startswith("THREAD-SECONDS"):
        return "thread_seconds"
    if title == "CONTINUOUS POSE WORK":
        return "cont_work"
    return "other"


def parse_file(path: Path, family: str, iteration: Optional[int], partition: Optional[int]) -> List[Record]:
    """Split a report file into records.

    A BENCHMARK RECORD header starts a record (the BENCH/ layout). In the files of
    earlier builds a file is one record, except that a repeated BENCHMARK CONTEXT
    header starts the next pass of the same file (the REC3D reports).
    """
    records: List[Record] = []
    current: Optional[Record] = None
    kind = ""
    saw_record_header = False

    def start(rec_family: str) -> Record:
        rec = Record(family=rec_family, source=path.name, iteration=iteration, partition=partition)
        records.append(rec)
        return rec

    for raw in path.read_text(errors="replace").splitlines():
        header = HEADER_RE.match(raw)
        if header:
            kind = section_kind(header.group(1))
            if kind == "record":
                saw_record_header = True
                current = start(family)
            elif kind == "context" and not saw_record_header and current is not None and \
                    (current.context or current.seconds):
                current = start(family)
            elif current is None:
                current = start(family)
            continue
        match = VALUE_RE.match(raw)
        if not match or current is None:
            continue
        key, text = " ".join(match.group(1).split()), match.group(2)
        value = as_number(text)
        if kind == "record":
            if key == "record family":
                current.family = text.strip()
            elif key == "record iteration" and value is not None:
                current.iteration = int(value)
            elif key == "record partition" and value is not None:
                current.partition = int(value)
        elif kind == "context":
            current.context[key] = text.strip()
        elif kind == "seconds" and value is not None:
            if "%" in key:
                current.percent[key] = value
            elif key.endswith("thread-s"):
                current.thread_seconds[key] = value
            else:
                current.seconds[key] = value
        elif kind == "percent" and value is not None:
            current.percent[key] = value
        elif kind == "thread_seconds" and value is not None:
            current.thread_seconds[key] = value
        elif kind == "cont_work" and value is not None:
            current.cont_work[key] = value
    records = [rec for rec in records if not rec.empty()]
    # passes of one family, iteration and partition, in file order
    seen: Dict[Tuple[str, Optional[int], Optional[int]], int] = defaultdict(int)
    for rec in records:
        seen[(rec.family, rec.iteration, rec.partition)] += 1
        rec.npass = seen[(rec.family, rec.iteration, rec.partition)]
    return records


def find_reports(indir: Path) -> Tuple[List[Path], str]:
    """The report files and their layout: 'BENCH/' or 'per-iteration files'."""
    bench = indir / "BENCH"
    if bench.is_dir() and any(bench.glob("*_BENCH*.txt")):
        return sorted(bench.glob("*_BENCH*.txt")), "BENCH/"
    files = sorted(indir.glob("*_BENCH*.txt"))
    if any(OLD_NAME_RE.match(p.name) and "_ITER" in p.name for p in files):
        return files, "per-iteration files"
    return files, "BENCH/"


def load_records(indir: Path) -> Tuple[List[Record], str, int]:
    files, layout = find_reports(indir)
    records: List[Record] = []
    for path in files:
        match = OLD_NAME_RE.match(path.name) if layout == "per-iteration files" else NEW_NAME_RE.match(path.name)
        if match is None:
            continue
        iteration = int(match.group("iter")) if "iter" in match.groupdict() and match.group("iter") else None
        partition = int(match.group("part")) if match.group("part") else None
        family = match.group("family")
        if layout == "per-iteration files" and family == "CAVGASSEMBLE" and iteration is None:
            family = "CAVGASSEMBLE_SUMS"   # the sum of the parts' contributions, named so since
        records.extend(parse_file(path, family, iteration, partition))
    return records, layout, len(files)


# ---------------------------------------------------------------- tables

def write_records_csv(records: Sequence[Record], outpath: Path) -> None:
    with outpath.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["family", "iteration", "partition", "pass", "section", "key", "value"])
        for rec in records:
            ident = [rec.family, "" if rec.iteration is None else rec.iteration,
                     "" if rec.partition is None else rec.partition, rec.npass]
            for section, values in (("context", rec.context), ("seconds", rec.seconds),
                                    ("percent", rec.percent), ("thread_seconds", rec.thread_seconds),
                                    ("cont_work", rec.cont_work)):
                for key, value in values.items():
                    writer.writerow(ident + [section, key, value])


def by_family(records: Sequence[Record], family: str) -> List[Record]:
    return [rec for rec in records if rec.family == family]


def first_pass_by_iteration(records: Sequence[Record]) -> Dict[int, Record]:
    out: Dict[int, Record] = {}
    for rec in records:
        if rec.iteration is not None and rec.npass == 1 and rec.iteration not in out:
            out[rec.iteration] = rec
    return out


def partition_stats(records: Sequence[Record], key: str) -> Dict[int, Tuple[float, float, float, int]]:
    """Per iteration: slowest, mean and fastest value of key over the partitions (the passes of a
    partition summed), and their number."""
    sums: Dict[Tuple[int, Optional[int]], float] = defaultdict(float)
    for rec in records:
        if rec.iteration is None or key not in rec.seconds:
            continue
        sums[(rec.iteration, rec.partition)] += rec.seconds[key]
    values: Dict[int, List[float]] = defaultdict(list)
    for (it, _), value in sums.items():
        values[it].append(value)
    return {it: (max(v), statistics.fmean(v), min(v), len(v)) for it, v in values.items() if v}


@dataclass
class Iteration:
    it: int
    mode: str = ""
    nspace: str = ""
    kto: str = ""
    scan: Optional[float] = None                                 # mean % of the search space scanned
    total: float = 0.0
    phases: Dict[str, float] = field(default_factory=dict)       # strategy phases as reported
    path: Dict[str, float] = field(default_factory=dict)         # critical-path split
    extra: Dict[str, float] = field(default_factory=dict)        # partition statistics and the rest


PATH_COMPONENTS = (
    ("setup", "setup/init", "#6f7d8c"),
    ("tables", "probability tables", "#2b6cb0"),
    ("prestep_overhead", "pre-step outside the tables", "#90cdf4"),
    ("search", "search", "#c05621"),
    ("sched_overhead", "scheduling outside the search", "#fbd38d"),
    ("volassemble", "volume assembly", "#2f855a"),
    ("assembly_overhead", "assembly outside volassemble", "#9ae6b4"),
    ("sigma2", "sigma2 commit", "#805ad5"),
    ("unaccounted", "unaccounted", "#e2e8f0"),
)


def build_iterations(records: Sequence[Record]) -> List[Iteration]:
    strategy = first_pass_by_iteration(by_family(records, "REFINE3D_STRATEGY"))
    if not strategy:
        return []
    volassemble = first_pass_by_iteration(by_family(records, "VOLASSEMBLE"))
    worker = partition_stats(by_family(records, "REFINE3D"), "match3D total time")
    sigma2_write = partition_stats(by_family(records, "REFINE3D"), "match3D sigma2 range write")
    tables: Dict[int, Tuple[float, float, float, int]] = {}
    fills: Dict[int, Tuple[float, float, float, int]] = {}
    for family in ("PROB_TAB", "PROB_TAB_NEIGH"):
        tables.update(partition_stats(by_family(records, family), "prob_tab total time"))
        fills.update(partition_stats(by_family(records, family), "prob_tab table fill"))
    rec3D = partition_stats(by_family(records, "REC3D"), "rec3D total time")
    rows: List[Iteration] = []
    for it in sorted(strategy):
        rec = strategy[it]
        row = Iteration(it=it, mode=rec.context.get("refine3D refine mode", ""),
                        nspace=rec.context.get("refine3D nspace", ""), kto=rec.context.get("refine3D kto", ""))
        row.scan = as_number(rec.context.get("refine3D % search space scanned", ""))
        s = rec.seconds
        row.total = s.get("refine3D total time", 0.0)
        row.phases = {
            "setup": s.get("refine3D strategy setup/init", 0.0),
            "model": s.get("refine3D reprojection model", 0.0),
            "prestep": s.get("refine3D probabilistic pre-step", 0.0),
            "sched": s.get("refine3D matcher/scheduler", 0.0),
            "assembly": s.get("refine3D assembly/postprocess", 0.0),
            "sigma2": s.get("refine3D sigma2 commit", 0.0),   # after the assembly, outside the four phases
        }
        tab_max = tables.get(it, (0.0, 0.0, 0.0, 0))[0]
        wrk_max = worker.get(it, (0.0, 0.0, 0.0, 0))[0]
        vol = volassemble[it].seconds.get("volassemble total time", 0.0) if it in volassemble else 0.0
        p = row.phases
        tab_on_path = min(tab_max, p["prestep"])
        wrk_on_path = min(wrk_max, p["sched"])
        vol_on_path = min(vol, p["assembly"])
        row.path = {
            "setup": p["setup"],
            "tables": tab_on_path,
            "prestep_overhead": p["prestep"] - tab_on_path,
            "search": wrk_on_path,
            "sched_overhead": p["sched"] - wrk_on_path,
            "volassemble": vol_on_path,
            "assembly_overhead": p["assembly"] - vol_on_path,
            "sigma2": p["sigma2"],
        }
        row.path["unaccounted"] = max(0.0, row.total - sum(row.path.values()))
        for name, stats in (("tables", tables), ("fill", fills), ("worker", worker), ("rec3D", rec3D),
                            ("sigma2_write", sigma2_write)):
            if it in stats:
                row.extra[f"{name}_max"], row.extra[f"{name}_mean"], row.extra[f"{name}_min"], \
                    row.extra[f"{name}_nparts"] = stats[it]
        if it in volassemble:
            row.extra["volassemble_nu_filter"] = volassemble[it].seconds.get("volassemble nonuniform_filter", 0.0)
            for key, label in NU_SPLIT:
                row.extra[f"volassemble_nu_{key}"] = volassemble[it].seconds.get(f"volassemble nu filter {label}", 0.0)
        rows.append(row)
    return rows


PATH_COMPONENTS_2D = (
    ("setup", "setup", "#6f7d8c"),
    ("tables", "2D tables", "#2b6cb0"),
    ("prestep_overhead", "pre-step outside the tables", "#90cdf4"),
    ("search", "search", "#c05621"),
    ("sched_overhead", "scheduling outside the search", "#fbd38d"),
    ("cavgassemble", "class-average assembly", "#2f855a"),
    ("assembly_overhead", "assembly outside cavgassemble", "#9ae6b4"),
    ("merge", "merge of the documents", "#d69e2e"),
    ("proj_io", "project writes and rereads", "#e53e3e"),
    ("sigma2", "sigma2 preparation and commit", "#805ad5"),
    ("finish", "convergence, star file, JPEG", "#b794f4"),
    ("unaccounted", "unaccounted", "#e2e8f0"),
)


def build_iterations_2D(records: Sequence[Record]) -> List[Iteration]:
    """One row per refine2D iteration with a REFINE2D_STRATEGY record: its phases and their critical path."""
    strategy = first_pass_by_iteration(by_family(records, "REFINE2D_STRATEGY"))
    if not strategy:
        return []
    search_recs = by_family(records, "REFINE2D")
    worker = partition_stats(search_recs, "match2D total time")
    tables = partition_stats(by_family(records, "PROB_TAB2D"), "prob_tab total time")
    fills = partition_stats(by_family(records, "PROB_TAB2D"), "prob_tab table fill")
    cavgasm = first_pass_by_iteration(by_family(records, "CAVGASSEMBLE"))
    first_search: Dict[int, Record] = {}
    for rec in search_recs:
        if rec.iteration is not None and rec.iteration not in first_search:
            first_search[rec.iteration] = rec
    rows: List[Iteration] = []
    for it in sorted(strategy):
        s = strategy[it].seconds
        ctx = first_search[it].context if it in first_search else {}
        row = Iteration(it=it, mode=ctx.get("match2D refine mode", ""), kto=ctx.get("match2D kto", ""))
        row.total = s.get("refine2D total time", 0.0)
        p = {key: s.get(f"refine2D {label}", 0.0) for key, label in (
            ("setup", "iteration setup"), ("prestep", "probabilistic pre-step"), ("sched", "matcher/scheduler"),
            ("merge", "merge alignment documents"), ("assembly", "class-average assembly"),
            ("proj_io", "project writes and rereads"), ("conv", "convergence"), ("star", "cls2D star file"),
            ("jpeg", "class-average JPEG"))}
        p["sigma2"] = s.get("refine2D sigma2 candidate prep", 0.0) + s.get("refine2D sigma2 commit", 0.0)
        row.phases = p
        tab_on_path = min(tables.get(it, (0.0, 0.0, 0.0, 0))[0], p["prestep"])
        wrk_on_path = min(worker.get(it, (0.0, 0.0, 0.0, 0))[0], p["sched"])
        cav = cavgasm[it].seconds.get("cavgassemble total time", 0.0) if it in cavgasm else 0.0
        cav_on_path = min(cav, p["assembly"])
        row.path = {
            "setup": p["setup"],
            "tables": tab_on_path,
            "prestep_overhead": p["prestep"] - tab_on_path,
            "search": wrk_on_path,
            "sched_overhead": p["sched"] - wrk_on_path,
            "cavgassemble": cav_on_path,
            "assembly_overhead": p["assembly"] - cav_on_path,
            "merge": p["merge"],
            "proj_io": p["proj_io"],
            "sigma2": p["sigma2"],
            "finish": p["conv"] + p["star"] + p["jpeg"],
        }
        row.path["unaccounted"] = max(0.0, row.total - sum(row.path.values()))
        for name, stats in (("tables", tables), ("fill", fills), ("worker", worker)):
            if it in stats:
                row.extra[f"{name}_max"], row.extra[f"{name}_mean"], row.extra[f"{name}_min"], \
                    row.extra[f"{name}_nparts"] = stats[it]
        rows.append(row)
    return rows


ITERATION_COLUMNS_2D = (
    ["iteration", "refine_mode", "kto", "strategy_total_s"]
    + [f"path_{k}_s" for k, _, _ in PATH_COMPONENTS_2D]
    + [f"{name}_{stat}_s" for name in ("tables", "fill", "worker") for stat in ("max", "mean", "min")]
    + ["nparts"]
)


def write_iterations_csv_2D(rows: Sequence[Iteration], outpath: Path) -> None:
    with outpath.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(ITERATION_COLUMNS_2D)
        for row in rows:
            values = [row.it, row.mode, row.kto, f"{row.total:.3f}"]
            values += [f"{row.path[k]:.3f}" for k, _, _ in PATH_COMPONENTS_2D]
            for name in ("tables", "fill", "worker"):
                for stat in ("max", "mean", "min"):
                    value = row.extra.get(f"{name}_{stat}")
                    values.append("" if value is None else f"{value:.3f}")
            nparts = int(row.extra.get("worker_nparts", row.extra.get("tables_nparts", 0)))
            values.append(nparts or "")
            writer.writerow(values)


def write_2D_path(rows: Sequence[Iteration], add) -> None:
    """The refine2D iterations' critical path, per stage (contiguous iterations of one refine mode) and the
    partitions' load balance."""
    if not rows:
        return
    total = sum(row.total for row in rows)
    add(f"refine2D: {len(rows)} iterations ({rows[0].it}-{rows[-1].it}); strategy wall time {total:.1f} s; "
        "its critical path:")
    add(f"  {'component':32s} {'seconds':>9s} {'share':>6s}")
    for key, label, _ in PATH_COMPONENTS_2D:
        value = sum(row.path[key] for row in rows)
        add(f"  {label:32s} {value:9.1f} {pct(value, total)}")
    add("")
    add("By stage (contiguous iterations with one refine mode), seconds per iteration:")
    short = (("setup", "setup"), ("tables", "tables"), ("prestep_overhead", "pre-ovh"), ("search", "search"),
             ("sched_overhead", "sch-ovh"), ("cavgassemble", "cavgasm"), ("assembly_overhead", "asm-ovh"),
             ("merge", "merge"), ("proj_io", "proj-io"), ("sigma2", "sigma2"), ("finish", "finish"))
    add(f"  {'iterations':>10s} {'mode':>10s} {'total':>7s} " + " ".join(f"{label:>7s}" for _, label in short))
    for block in stages_of(rows):
        n = len(block)
        add(f"  {block[0].it:>4d}-{block[-1].it:<5d} {block[0].mode:>10s} {sum(r.total for r in block) / n:7.2f} "
            + " ".join(f"{sum(r.path[k] for r in block) / n:7.2f}" for k, _ in short))
    add("")
    add("refine2D load balance (slowest partition over the mean, median over iterations):")
    found = False
    for name, label in (("tables", "2D table parts"), ("worker", "search partitions")):
        ratios = [row.extra[f"{name}_max"] / row.extra[f"{name}_mean"] for row in rows
                  if row.extra.get(f"{name}_mean", 0.0) > 0 and row.extra.get(f"{name}_nparts", 0) > 1]
        if ratios:
            found = True
            add(f"  {label:26s} median {statistics.median(ratios):.2f}, max {max(ratios):.2f} "
                f"({len(ratios)} iterations)")
    if not found:
        add("  one partition")
    add("")


# the nonuniform filter's split in the volassemble report (column suffix, report label)
NU_SPLIT = (("candidate_filter", "candidate filter"), ("objective", "objective"), ("optimisation", "optimisation"),
            ("composition", "composition"))

ITERATION_COLUMNS = (
    ["iteration", "refine_mode", "nspace", "kto", "strategy_total_s"]
    + [f"strategy_{k}_s" for k in ("setup", "model", "prestep", "sched", "assembly", "sigma2")]
    + [f"path_{k}_s" for k, _, _ in PATH_COMPONENTS]
    + [f"{name}_{stat}_s" for name in ("tables", "fill", "worker", "rec3D", "sigma2_write")
       for stat in ("max", "mean", "min")]
    + ["nparts", "volassemble_nu_filter_s", "search_space_scanned_pct"]
    + [f"volassemble_nu_{key}_s" for key, _ in NU_SPLIT]
)


def write_iterations_csv(rows: Sequence[Iteration], outpath: Path) -> None:
    with outpath.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(ITERATION_COLUMNS)
        for row in rows:
            values = [row.it, row.mode, row.nspace, row.kto, f"{row.total:.3f}"]
            values += [f"{row.phases[k]:.3f}" for k in ("setup", "model", "prestep", "sched", "assembly", "sigma2")]
            values += [f"{row.path[k]:.3f}" for k, _, _ in PATH_COMPONENTS]
            for name in ("tables", "fill", "worker", "rec3D", "sigma2_write"):
                for stat in ("max", "mean", "min"):
                    value = row.extra.get(f"{name}_{stat}")
                    values.append("" if value is None else f"{value:.3f}")
            nparts = int(row.extra.get("worker_nparts", row.extra.get("tables_nparts", 0)))
            values += [nparts or "", f"{row.extra.get('volassemble_nu_filter', 0.0):.3f}",
                       "" if row.scan is None else f"{row.scan:.3f}"]
            values += [f"{row.extra.get(f'volassemble_nu_{key}', 0.0):.3f}" for key, _ in NU_SPLIT]
            writer.writerow(values)


# ---------------------------------------------------------------- report

def stages_of(rows: Sequence[Iteration]) -> List[List[Iteration]]:
    """Contiguous iterations with the same refine mode and nspace."""
    blocks: List[List[Iteration]] = []
    for row in rows:
        if blocks and (blocks[-1][-1].mode, blocks[-1][-1].nspace) == (row.mode, row.nspace):
            blocks[-1].append(row)
        else:
            blocks.append([row])
    return blocks


def mad(values: Sequence[float]) -> float:
    med = statistics.median(values)
    return statistics.median([abs(v - med) for v in values])


def pct(part: float, whole: float) -> str:
    return f"{100.0 * part / whole:5.1f}%" if whole > 0 else "    -"


def phase_means(records: Sequence[Record]) -> List[Tuple[str, float, float]]:
    """Mean and maximum of every timing over the records, largest mean first."""
    values: Dict[str, List[float]] = defaultdict(list)
    for rec in records:
        for key, value in rec.seconds.items():
            values[key].append(value)
    out = [(key, statistics.fmean(v), max(v)) for key, v in values.items()]
    return sorted(out, key=lambda item: -item[1])


TABLE_FAMILIES = (("PROB_TAB", "prob"), ("PROB_TAB_NEIGH", "prob_neigh"), ("PROB_TAB2D", "2D"))
THREAD_FAMILIES = TABLE_FAMILIES + (("REFINE2D", "2D search and class sums"), ("REFINE3D", "3D search"),
                                    ("CAVGSUMS", "offline class sums"))


def write_thread_seconds(records: Sequence[Record], add) -> None:
    """The THREAD-SECONDS blocks of the table parts, the 2D search and the class sums, summed over parts and
    iterations per family: the additive entries summed, the ratios (load balance, cost per score) and the
    team sizes averaged over the records."""
    for family, label in THREAD_FAMILIES:
        recs = [rec for rec in by_family(records, family) if rec.thread_seconds]
        if not recs:
            continue
        sums: Dict[str, float] = defaultdict(float)
        ratios: Dict[str, List[float]] = defaultdict(list)
        for rec in recs:
            for key, value in rec.thread_seconds.items():
                if "load balance" in key or "/(" in key or key.endswith("threads"):
                    ratios[key].append(value)
                else:
                    sums[key] += value
        add(f"Thread-seconds ({label}), summed over {len(recs)} records:")
        for key, value in sums.items():
            add(f"  {key:46s} {value:10.1f}")
        for key, values in ratios.items():
            add(f"  {key:46s} {statistics.fmean(values):10.4g} (mean)")
        add("")


CONT_PREFIX = "match3D cont "
CONT_COUNTS = ("solves", "accepted", "no improvement", "no reliable update", "bound rejected", "invalid numerics",
               "iteration limit", "invalid preparation", "not attempted", "solves at the iteration cap",
               "LM iterations", "trial evaluations", "accepted steps", "bound hits", "calculator passes",
               "pixel-evaluations")
CONT_RATIO_BINS = ("0.25-0.75", "0.75-1.25", "1.25-1.75", "1.75-2.5", ">= 2.5")


def write_cont_work(records: Sequence[Record], add) -> None:
    """The CONTINUOUS POSE WORK sections of the 3D search records (refine=cont and the polish), summed over
    parts and iterations; the per-solve figures and the histogram of actual over predicted reduction at the
    accepted steps are recomputed from the sums, the time per pixel-evaluation from the search's
    thread-seconds."""
    recs = [rec for rec in by_family(records, "REFINE3D") if rec.cont_work]
    if not recs:
        return
    sums = {name: sum(rec.cont_work.get(CONT_PREFIX + name, 0.0) for rec in recs) for name in CONT_COUNTS}
    bins = {name: sum(rec.cont_work.get(CONT_PREFIX + "actual/predicted " + name, 0.0) for rec in recs)
            for name in CONT_RATIO_BINS}
    search_ts = sum(rec.thread_seconds.get("match3D search", 0.0) for rec in recs)
    # a solve that evaluated its seed; the two statuses without one make no pass
    solved = sums["solves"] - sums["not attempted"] - sums["invalid preparation"]
    add(f"Continuous pose work (refine=cont and the polish), summed over {len(recs)} records:")
    for name in CONT_COUNTS:
        add(f"  {name:46s} {sums[name]:14.0f}")
    if solved > 0:
        add(f"  {'passes per solve':46s} {sums['calculator passes'] / solved:14.3f}")
        add(f"  {'LM iterations per solve':46s} {sums['LM iterations'] / solved:14.3f}")
    if sums["trial evaluations"] > 0:
        add(f"  {'rejected fraction of trial evaluations':46s} "
            f"{1.0 - sums['accepted steps'] / sums['trial evaluations']:14.3f}")
    nbinned = sum(bins.values())
    if nbinned > 0:
        for name in CONT_RATIO_BINS:
            add(f"  {'accepted steps with actual/predicted ' + name:46s} {bins[name] / nbinned:14.3f}")
    if sums["pixel-evaluations"] > 0 and search_ts > 0:
        add(f"  {'search ns per pixel-evaluation':46s} {1e9 * search_ts / sums['pixel-evaluations']:14.2f}")
    add("")


PRESTEP_KEYS = (("scheduling (parts)", "parts"), ("table reads", "reads"), ("ref_assign", "assign"),
                ("total time", "total"))


def write_prestep_groups(records: Sequence[Record], add) -> None:
    """The probabilistic pre-step driver (PROB_ALIGN, PROB_ALIGN_NEIGH) per table, number of states and
    directions: the mean seconds per iteration of the parts' scheduling, the master's table reads, its
    ref_assign and the total."""
    groups: Dict[Tuple[str, str, str], List[Record]] = defaultdict(list)
    for family, label in (("PROB_ALIGN", "prob"), ("PROB_ALIGN_NEIGH", "prob_neigh")):
        for rec in first_pass_by_iteration(by_family(records, family)).values():
            groups[(label, rec.context.get("prob_align states", "?"),
                    rec.context.get("prob_align projection directions", "?"))].append(rec)
    if not groups:
        return
    add("Probabilistic pre-step driver by table, states and directions, mean seconds per iteration:")
    add(f"  {'table':>10s} {'states':>6s} {'nspace':>6s} {'iters':>5s} "
        + " ".join(f"{short:>7s}" for _, short in PRESTEP_KEYS))
    for (label, nstates, nspace), recs in groups.items():
        means = [statistics.fmean(rec.seconds.get(f"prob_align {key}", 0.0) for rec in recs)
                 for key, _ in PRESTEP_KEYS]
        add(f"  {label:>10s} {nstates:>6s} {nspace:>6s} {len(recs):5d} " + " ".join(f"{m:7.3f}" for m in means))
    add("")


def write_2D_tables(records: Sequence[Record], add) -> None:
    """The refine2D iterations with a 2D probability table: the slowest table part per iteration (total and
    fill) summed, against the refine2D strategy's wall time of the same iterations."""
    strategy = first_pass_by_iteration(by_family(records, "REFINE2D_STRATEGY"))
    tables = partition_stats(by_family(records, "PROB_TAB2D"), "prob_tab total time")
    fills = partition_stats(by_family(records, "PROB_TAB2D"), "prob_tab table fill")
    its = sorted(it for it in tables if it in strategy)
    if not its:
        return
    total = sum(strategy[it].seconds.get("refine2D total time", 0.0) for it in its)
    slowest = sum(tables[it][0] for it in its)
    fill = sum(fills[it][0] for it in its if it in fills)
    ratios = [tables[it][0] / tables[it][1] for it in its if tables[it][1] > 0]
    add(f"refine2D iterations with a 2D table: {len(its)}; their refine2D wall time {total:.1f} s")
    add(f"  slowest table part, total                  {slowest:9.1f} s  {pct(slowest, total)}")
    add(f"  slowest table part, table fill             {fill:9.1f} s  {pct(fill, total)}")
    if ratios:
        add(f"  table parts' load balance (slowest over the mean): median {statistics.median(ratios):.2f}")
    add("")


def write_report(records: Sequence[Record], rows: Sequence[Iteration], rows2D: Sequence[Iteration], layout: str,
                 nfiles: int, indir: Path, outpath: Path) -> None:
    lines: List[str] = []
    add = lines.append
    add("SIMPLE benchmark summary")
    add("=" * 24)
    add(f"input   : {indir}")
    add(f"layout  : {layout}, {nfiles} files, {len(records)} records")
    counts: Dict[str, int] = defaultdict(int)
    keys: Dict[Tuple[str, Optional[int], Optional[int]], int] = defaultdict(int)
    for rec in records:
        counts[rec.family] += 1
        if rec.family != "REC3D" and rec.iteration is not None:
            keys[(rec.family, rec.iteration, rec.partition)] += 1
    duplicates = sum(1 for n in keys.values() if n > 1)
    add("families: " + ", ".join(f"{fam} {n}" for fam, n in sorted(counts.items())))
    if duplicates:
        add(f"WARNING : {duplicates} family/iteration/partition keys occur more than once, as after a second "
            "run in the same directory; the first record of each is summarized")
    add("")

    for rec in by_family(records, "SOLVE3D"):
        c, s = rec.context, rec.seconds
        total = s.get("solve3D total time", 0.0)
        add(f"Run: {c.get('solve3D route', 'solve3D').strip()}, nstates {c.get('solve3D nstates', '?')}, "
            f"box {c.get('solve3D box', '?')}, stages {c.get('solve3D first stage', '?')}-"
            f"{c.get('solve3D last stage', '?')}, {c.get('solve3D searched particles', '?')} particles searched, "
            f"nsample {c.get('solve3D nsample', '?')}")
        add(f"  total {total:.1f} s")
        for key, value in s.items():
            if key == "solve3D total time":
                continue
            add(f"  {key.replace('solve3D ', ''):34s} {value:9.1f} s  {pct(value, total)}")
        add("")

    if not rows:
        add("No refine3D strategy records: the critical path and the figure need REFINE3D_STRATEGY reports.")
        add("")
    else:
        its = [row.it for row in rows]
        missing = sorted(set(range(its[0], its[-1] + 1)) - set(its))
        total = sum(row.total for row in rows)
        add(f"refine3D: {len(rows)} iterations ({its[0]}-{its[-1]}"
            + (f"; no strategy record for {', '.join(map(str, missing))}" if missing else "") + ")")
        add(f"strategy wall time {total:.1f} s; its critical path:")
        add(f"  {'component':32s} {'seconds':>9s} {'share':>6s}")
        for key, label, _ in PATH_COMPONENTS:
            value = sum(row.path[key] for row in rows)
            add(f"  {label:32s} {value:9.1f} {pct(value, total)}")
        compute = sum(row.path["tables"] + row.path["search"] + row.path["volassemble"] for row in rows)
        overhead = sum(row.path["prestep_overhead"] + row.path["sched_overhead"] + row.path["assembly_overhead"]
                       for row in rows)
        add(f"  slowest-partition and assembly compute {compute:.1f} s ({pct(compute, total).strip()}); "
            f"master time around it {overhead:.1f} s ({pct(overhead, total).strip()})")
        add("")
        add("By stage (contiguous iterations with one refine mode and nspace), seconds per iteration; scan% is the")
        add("mean percentage of the search space scanned at the stage's last iteration:")
        header = (f"  {'iterations':>10s} {'mode':>10s} {'nspace':>6s} {'total':>7s} {'setup':>6s} {'tables':>7s} "
                  f"{'pre-ovh':>7s} {'search':>7s} {'sch-ovh':>7s} {'volasm':>6s} {'asm-ovh':>7s} {'sigma2':>6s} "
                  f"{'scan%':>6s}")
        add(header)
        for block in stages_of(rows):
            n = len(block)

            def mean_of(key: str) -> float:
                return sum(r.path[key] for r in block) / n
            add(f"  {block[0].it:>4d}-{block[-1].it:<5d} {block[0].mode:>10s} {block[0].nspace:>6s} "
                f"{sum(r.total for r in block) / n:7.2f} {mean_of('setup'):6.2f} {mean_of('tables'):7.2f} "
                f"{mean_of('prestep_overhead'):7.2f} {mean_of('search'):7.2f} {mean_of('sched_overhead'):7.2f} "
                f"{mean_of('volassemble'):6.2f} {mean_of('assembly_overhead'):7.2f} {mean_of('sigma2'):6.2f} "
                + ("     -" if block[-1].scan is None else f"{block[-1].scan:6.2f}"))
        add("")
        add("Load balance (slowest partition over the mean, median over iterations):")
        if not any(row.extra.get("worker_mean", 0.0) > 0 or row.extra.get("tables_mean", 0.0) > 0 for row in rows):
            add("  no partition reports")
        for name, label in (("tables", "probability-table parts"), ("worker", "search workers"),
                            ("rec3D", "reconstruction passes")):
            ratios = [row.extra[f"{name}_max"] / row.extra[f"{name}_mean"] for row in rows
                      if row.extra.get(f"{name}_mean", 0.0) > 0]
            if ratios:
                add(f"  {label:26s} median {statistics.median(ratios):.2f}, max {max(ratios):.2f} "
                    f"({len(ratios)} iterations)")
        add("")
        add("Outlier iterations (total above the median of the iterations with the same refine mode, nspace and")
        add("low-pass shell by more than 3 MAD and 10%):")
        found = False
        groups: Dict[Tuple[str, str, str], List[Iteration]] = defaultdict(list)
        for row in rows:
            groups[(row.mode, row.nspace, row.kto)].append(row)
        for block in groups.values():
            if len(block) < 5:
                continue
            totals = [r.total for r in block]
            med, spread = statistics.median(totals), mad(totals)
            for row in block:
                if row.total > med + max(3.0 * spread, 0.10 * med):
                    medians = {k: statistics.median([r.path[k] for r in block]) for k, _, _ in PATH_COMPONENTS}
                    key = max(medians, key=lambda k: row.path[k] - medians[k])
                    label = dict((k, l) for k, l, _ in PATH_COMPONENTS)[key]
                    add(f"  iteration {row.it:3d}: {row.total:6.2f} s against a median of {med:.2f} s; "
                        f"mostly {label} ({row.path[key]:.2f} s against {medians[key]:.2f} s)")
                    found = True
        if not found:
            add("  none")
        add("")
    write_2D_path(rows2D, add)
    write_2D_tables(records, add)
    write_thread_seconds(records, add)
    write_cont_work(records, add)
    write_prestep_groups(records, add)
    for family, title in (("REFINE3D", "Search worker phases (mean and maximum over iterations and partitions)"),
                          ("PROB_TAB", "Probability-table part phases (prob)"),
                          ("PROB_TAB_NEIGH", "Probability-table part phases (prob_neigh)"),
                          ("PROB_ALIGN", "Probabilistic pre-step driver phases (prob)"),
                          ("PROB_ALIGN_NEIGH", "Probabilistic pre-step driver phases (prob_neigh)"),
                          ("REC3D", "Reconstruction pass phases"),
                          ("VOLASSEMBLE", "Volume assembly phases"),
                          ("SIGMA2_COMMIT", "Sigma2 commit phases"),
                          ("REFINE2D_STRATEGY", "refine2D strategy phases"),
                          ("REFINE2D", "2D search phases (every partition)"),
                          ("PROB_ALIGN2D", "2D probabilistic pre-step driver phases"),
                          ("SOLVE2D", "solve2D phases"),
                          ("PROB_TAB2D", "2D probability-table part phases"),
                          ("CAVGASSEMBLE", "Class-average assembly phases"),
                          ("CAVGASSEMBLE_SUMS", "Class-average sums of the parts"),
                          ("CAVGSUMS", "Offline class sums per partition")):
        recs = by_family(records, family)
        if not recs:
            continue
        add(f"{title}, {len(recs)} records:")
        for key, mean, maximum in phase_means(recs):
            add(f"  {key:46s} mean {mean:8.3f} s   max {maximum:8.3f} s")
        add("")
    outpath.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------- figure

def linear(value: float, low: float, high: float, pix_low: float, pix_high: float) -> float:
    if high <= low:
        return (pix_low + pix_high) / 2.0
    return pix_low + (value - low) * (pix_high - pix_low) / (high - low)


def nice_limit(values: Iterable[float]) -> float:
    maximum = max(values, default=1.0)
    if maximum <= 0:
        return 1.0
    step = 10 ** math.floor(math.log10(maximum))
    return math.ceil(maximum / step * 5.0) / 5.0 * step


def tick_step(limit: float) -> float:
    raw = limit / 5.0
    power = 10 ** math.floor(math.log10(raw))
    for mult in (1, 2, 2.5, 5, 10):
        if raw <= mult * power:
            return mult * power
    return 10 * power


def svg_text(x: float, y: float, text: str, size: int = 16, anchor: str = "start", weight: int = 400) -> str:
    return (f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-size="{size}" '
            f'font-weight="{weight}">{escape(text)}</text>')


def svg_vertical_text(x: float, y: float, text: str, size: int = 16) -> str:
    return (f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" font-size="{size}" '
            f'transform="rotate(-90 {x:.1f} {y:.1f})">{escape(text)}</text>')


def legend(items: Sequence[Tuple[str, str, str]], x: float, y: float, width: float) -> Tuple[List[str], float]:
    """Swatches ('box') or line samples ('line') wrapped to width; returns the elements and the next y."""
    out: List[str] = []
    xpos, ypos = x, y
    for label, color, style in items:
        item_w = (40 if style == "line" else 28) + len(label) * 10.5
        if xpos + item_w > x + width and xpos > x:
            xpos, ypos = x, ypos + 30
        if style == "line":
            out.append(f'<line x1="{xpos:.1f}" y1="{ypos - 6:.1f}" x2="{xpos + 28:.1f}" y2="{ypos - 6:.1f}" '
                       f'stroke="{color}" stroke-width="3"/>')
            out.append(svg_text(xpos + 36, ypos, label, 18))
        else:
            out.append(f'<rect x="{xpos:.1f}" y="{ypos - 15:.1f}" width="18" height="18" fill="{color}" '
                       f'stroke="#a0aec0" stroke-width="0.5"/>')
            out.append(svg_text(xpos + 26, ypos, label, 18))
        xpos += item_w
    return out, ypos + 30


def axes(out: List[str], plot_l: float, plot_r: float, top: float, bottom: float, limit: float,
         first: int, last: int, rows: Sequence[Iteration], ylabel: str) -> None:
    step = tick_step(limit)
    tick = 0.0
    while tick <= limit + 1e-9:
        y = linear(tick, 0, limit, bottom, top)
        out.append(f'<line class="grid" x1="{plot_l}" y1="{y:.1f}" x2="{plot_r}" y2="{y:.1f}"/>')
        out.append(svg_text(plot_l - 12, y + 6, f"{tick:g}", 18, "end"))
        tick += step
    out.append(f'<line class="axis" x1="{plot_l}" y1="{top}" x2="{plot_l}" y2="{bottom}"/>')
    out.append(f'<line class="axis" x1="{plot_l}" y1="{bottom}" x2="{plot_r}" y2="{bottom}"/>')
    out.append(svg_vertical_text(plot_l - 70, (top + bottom) / 2, ylabel, 20))
    for it in range(first, last + 1):
        if (it - first) % 10 != 0:
            continue
        x = linear(it, first, last, plot_l, plot_r)
        out.append(f'<line class="axis" x1="{x:.1f}" y1="{bottom}" x2="{x:.1f}" y2="{bottom + 6}"/>')
        out.append(svg_text(x, bottom + 28, str(it), 18, "middle"))
    for block in stages_of(rows)[1:]:
        x = linear(block[0].it - 0.5, first, last, plot_l, plot_r)
        out.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{bottom}" stroke="#1a202c" '
                   f'stroke-width="1.5" stroke-dasharray="5 4"/>')


def write_figure(rows: Sequence[Iteration], outpath: Path, title: str,
                 components: Sequence[Tuple[str, str, str]] = PATH_COMPONENTS, workflow: str = "refine3D") -> None:
    first, last = rows[0].it, rows[-1].it
    width, height = 1800, 1500
    plot_l, plot_r = 130, width - 60
    top_a, bottom_a = 200, 760
    top_b, bottom_b = 1070, 1380
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<title>{workflow} benchmark profile</title>',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#17202a}'
        '.axis{stroke:#4a5568;stroke-width:2}.grid{stroke:#e2e8f0;stroke-width:1}</style>',
        svg_text(plot_l, 55, title, 34, weight=600),
        svg_text(plot_l, 92, "A: the strategy's wall time per iteration, split into the slowest partition's work and "
                 "the master's time around it. B: the slowest and the mean partition.", 20),
    ]
    # panel A: stacked critical path
    limit_a = nice_limit(row.total for row in rows)
    out.append(svg_text(plot_l, top_a - 50, "A  Critical path of each iteration", 24, weight=600))
    for block in stages_of(rows):
        x = linear(block[0].it - 0.5, first, last, plot_l, plot_r)
        label = f"{block[0].mode} {block[0].nspace}"
        if x + 6 + 9 * len(label) > plot_r:
            out.append(svg_text(plot_r, top_a - 14, label, 16, "end"))
        else:
            out.append(svg_text(x + 6, top_a - 14, label, 16))
    axes(out, plot_l, plot_r, top_a, bottom_a, limit_a, first, last, rows, "wall time (s)")
    bar_w = max(2.0, (plot_r - plot_l) / max(1, last - first + 1) * 0.8)
    for row in rows:
        x = linear(row.it, first, last, plot_l, plot_r)
        y = bottom_a
        for key, _, color in components:
            h = bottom_a - linear(row.path[key], 0, limit_a, bottom_a, top_a)
            if h <= 0:
                continue
            y -= h
            out.append(f'<rect x="{x - bar_w / 2:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" '
                       f'fill="{color}"/>')
    out.append(svg_text((plot_l + plot_r) / 2, bottom_a + 58, "refinement iteration", 20, "middle"))
    items, next_y = legend([(label, color, "box") for _, label, color in components],
                           plot_l, bottom_a + 95, plot_r - plot_l)
    out += items
    # panel B: slowest and mean partition
    series = (("tables_max", "probability tables, slowest part", "#2b6cb0", ""),
              ("tables_mean", "probability tables, mean part", "#2b6cb0", "6 5"),
              ("worker_max", "search, slowest worker", "#c05621", ""),
              ("worker_mean", "search, mean worker", "#c05621", "6 5"))
    present = [s for s in series if any(s[0] in row.extra for row in rows)]
    if present:
        limit_b = nice_limit(row.extra.get(key, 0.0) for row in rows for key, _, _, _ in present)
        out.append(svg_text(plot_l, top_b - 50, "B  Partition wall time: slowest against mean", 24, weight=600))
        axes(out, plot_l, plot_r, top_b, bottom_b, limit_b, first, last, rows, "wall time per part (s)")
        for key, _, color, dash in present:
            points = [f"{linear(row.it, first, last, plot_l, plot_r):.1f},"
                      f"{linear(row.extra[key], 0, limit_b, bottom_b, top_b):.1f}" for row in rows if key in row.extra]
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
            out.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.5"{dash_attr} '
                       f'points="{" ".join(points)}"/>')
        out.append(svg_text((plot_l + plot_r) / 2, bottom_b + 58, "refinement iteration", 20, "middle"))
        items, _ = legend([(label, color, "line") for _, label, color, _ in present], plot_l, top_b - 18,
                          plot_r - plot_l)
        out += items
    out.append('</svg>')
    outpath.write_text("\n".join(out) + "\n")


def main(argv: List[str]) -> int:
    if not 2 <= len(argv) <= 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    indir = Path(argv[1]).expanduser().resolve()
    if not indir.is_dir():
        print(f"Benchmark directory does not exist: {indir}", file=sys.stderr)
        return 2
    records, layout, nfiles = load_records(indir)
    if not records:
        print(f"No benchmark reports in {indir}", file=sys.stderr)
        return 1
    reports_dir = indir / "BENCH" if (indir / "BENCH").is_dir() else indir
    outdir = Path(argv[2]).expanduser().resolve() if len(argv) == 3 else reports_dir / "bench_summary"
    outdir.mkdir(parents=True, exist_ok=True)
    rows = build_iterations(records)
    rows2D = build_iterations_2D(records)
    write_records_csv(records, outdir / "bench_records.csv")
    write_iterations_csv(rows, outdir / "bench_iterations.csv")
    write_report(records, rows, rows2D, layout, nfiles, indir, outdir / "bench_report.txt")
    written = ["bench_records.csv", "bench_iterations.csv", "bench_report.txt"]
    if rows:
        write_figure(rows, outdir / "bench_profile.svg", f"refine3D benchmark profile: {indir.name}")
        written.append("bench_profile.svg")
    if rows2D:
        write_iterations_csv_2D(rows2D, outdir / "bench_iterations2D.csv")
        write_figure(rows2D, outdir / "bench_profile2D.svg", f"refine2D benchmark profile: {indir.name}",
                     PATH_COMPONENTS_2D, "refine2D")
        written += ["bench_iterations2D.csv", "bench_profile2D.svg"]
    print(f"Parsed {len(records)} records from {nfiles} files ({layout}) in {indir}")
    for name in written:
        print(f"Wrote {outdir / name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
