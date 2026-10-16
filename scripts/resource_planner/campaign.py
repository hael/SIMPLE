"""Bounded local campaigns, immutable input identities, and append-only evidence."""

from contextlib import contextmanager
import csv
import hashlib
import json
import os
import platform
from pathlib import Path
import random
import shutil
import signal
import subprocess
import threading
import time

from .adapters import build_command, positive_int, validate_case


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def fingerprint(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_files(root):
    result = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"fixture symlinks are not supported: {path}")
        if path.is_file():
            result[str(path.relative_to(root))] = file_hash(path)
        elif not path.is_dir():
            raise ValueError(f"fixture contains a special file: {path}")
    return result


def relative_file(value):
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or path == Path("."):
        raise ValueError(f"expected a relative file within the run directory: {value}")
    return path


def argv(value):
    if not isinstance(value, list) or not value or any(
            not isinstance(item, str) or not item or "\x00" in item for item in value):
        raise ValueError("commands must be nonempty arrays of argument strings")
    return value


def load_plan(path):
    path = path.resolve()
    manifest = json.loads(path.read_text())
    # Reject nonfinite JSON anywhere, including opaque fixture features.
    canonical(manifest)
    if manifest["schema_version"] != 1 or manifest["backend"] != "local":
        raise ValueError("only campaign schema 1 and backend=local are supported")
    budget = manifest["budget"]
    for key in ("cpus", "wall_seconds", "cpu_seconds", "attempt_seconds", "disk_bytes"):
        positive_int(budget[key], key)
    available = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    if available is not None and budget["cpus"] > available:
        raise ValueError("CPU budget exceeds process affinity")
    # A seed controls scheduling and is propagated independently of repeat order.
    positive_int(manifest["seed"], "seed")
    if positive_int(manifest["repetitions"], "repetitions") < 3:
        raise ValueError("at least three measured repetitions are required")
    fixture = (path.parent / manifest["fixture"]).resolve()
    if not fixture.is_dir():
        raise ValueError("fixture must be an existing directory")
    files = tree_files(fixture)
    if not files:
        raise ValueError("fixture must not be empty")
    fixture_bytes = sum((fixture / name).stat().st_size for name in files)
    executables = {}
    for key, value in manifest["executables"].items():
        executable = (path.parent / value).resolve()
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise ValueError(f"executable unavailable: {executable}")
        executables[key] = {"path": str(executable), "sha256": file_hash(executable)}
    cases, ids, runs = manifest["cases"], set(), []
    if not isinstance(cases, list) or not cases:
        raise ValueError("cases must be nonempty")
    for index, case in enumerate(cases):
        if not isinstance(case["id"], str) or not case["id"] or case["id"] in ids:
            raise ValueError("case IDs must be nonempty and unique")
        ids.add(case["id"])
        validate_case(case, budget["cpus"])
        relative_file(case["project"])
        if not isinstance(case["prepare"], list) or not case["prepare"]:
            raise ValueError("each case requires fixture preparation commands")
        for command in case["prepare"]:
            argv(command)
        argv(case["validate"])
        for command in case["prepare"] + [case["validate"]]:
            tool = shutil.which(command[0])
            if not tool:
                raise ValueError(f"hook executable unavailable: {command[0]}")
            executables["hook:" + command[0]] = {"path": str(Path(tool).resolve()),
                                                 "sha256": file_hash(Path(tool))}
        key = "private" if case["mode"] == "worker" else "simple"
        if key not in executables:
            raise ValueError(f"missing {key} executable")
        seen = set()
        for candidate in case["candidates"]:
            identity = fingerprint(candidate)
            if identity in seen:
                raise ValueError("duplicate candidate")
            seen.add(identity)
        # Warmup each candidate, then randomize the measured runs in blocks.
        for repeat in range(manifest["repetitions"] + 1):
            candidates = list(case["candidates"])
            random.Random(manifest["seed"] + index * 1000 + repeat).shuffle(candidates)
            for candidate in candidates:
                runs.append({"case_index": index, "candidate": candidate,
                             "repeat": repeat, "warmup": repeat == 0})
    seconds = len(runs) * budget["attempt_seconds"]
    if seconds > budget["wall_seconds"] or seconds * budget["cpus"] > budget["cpu_seconds"]:
        raise ValueError("worst-case campaign exceeds wall/CPU reservation budget; reduce candidates or timeout")
    if fixture_bytes * len(runs) > budget["disk_bytes"]:
        raise ValueError("fixture copies alone exceed campaign disk budget")
    cpu_ids = sorted(os.sched_getaffinity(0))[:budget["cpus"]] if hasattr(os, "sched_getaffinity") else None
    environment = {key: value for key, value in os.environ.items()
                   if key.startswith(("OMP_", "MKL_", "OPENBLAS_", "SIMPLE_")) or
                   key in ("PATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH")}
    plan = {"manifest": manifest, "fixture_files": files, "fixture": str(fixture),
            "executables": executables, "runs": runs, "cpu_ids": cpu_ids,
            "environment": environment, "platform": platform.platform(),
            "implementation": {name: file_hash(Path(__file__).with_name(name))
                               for name in ("campaign.py", "adapters.py")}}
    plan["campaign_id"] = fingerprint(plan)
    for run in runs:
        run["experiment_id"] = fingerprint({"campaign": plan["campaign_id"], **run})
    return plan


def render(plan, destination):
    output = []
    for run in plan["runs"]:
        case = plan["manifest"]["cases"][run["case_index"]]
        directory = destination / run["experiment_id"] / "attempt-1"
        key = "private" if case["mode"] == "worker" else "simple"
        output.append({**run, "case": case["id"], "command": build_command(
            case, run["candidate"], plan["executables"][key]["path"], directory, directory / "fixture")})
    return {"campaign_id": plan["campaign_id"], "budget": plan["manifest"]["budget"],
            "execution_enabled": False, "runs": output}


def append_record(path, record):
    with path.open("a", encoding="utf-8") as stream:
        stream.write(canonical(record).decode() + "\n")
        stream.flush()
        os.fsync(stream.fileno())


@contextmanager
def campaign_lock(destination):
    # Advisory lock releases on crashes; no stale PID lock to remove on resume.
    import fcntl
    with (destination / "campaign.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("another runner owns this campaign") from error
        yield


def disk_usage(directory):
    return sum(path.stat().st_size for path in directory.rglob("*")
               if path.is_file() and not path.is_symlink())


def run_process(command, cwd, env, deadline, stdout, stderr, destination, disk_limit, cpu_ids=None):
    started = time.monotonic()
    with stdout.open("wb") as out, stderr.open("wb") as err:
        bind = (lambda: os.sched_setaffinity(0, cpu_ids)) if cpu_ids is not None else None
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=out, stderr=err,
                                   start_new_session=True, preexec_fn=bind)
        stopped = threading.Event()
        reasons = []
        def kill_group():
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        def monitor():
            try:
                while not stopped.is_set():
                    if time.monotonic() >= deadline:
                        reasons.append("timeout")
                        kill_group()
                        return
                    if disk_usage(destination) > disk_limit:
                        reasons.append("disk_budget_exceeded")
                        kill_group()
                        return
                    stopped.wait(0.1)
            except OSError as error:
                reasons.append(f"resource_monitor_failed: {error}")
                kill_group()
        watcher = threading.Thread(target=monitor, daemon=True)
        watcher.start()
        try:
            process.wait()
            # Blocking wait avoids adding the resource polling interval to timings.
            elapsed = time.monotonic() - started
        finally:
            stopped.set()
            # Also remove descendants left behind after a launcher exits.
            kill_group()
            process.wait()
            watcher.join()
    reason = reasons[0] if reasons else None
    if disk_usage(destination) > disk_limit:
        reason = "disk_budget_exceeded"
    return {"returncode": process.returncode, "reason": reason,
            "wall_seconds": elapsed}


def telemetry(directory):
    peaks, errors = [], []
    for path in sorted(directory.glob("memory_usage_*.csv")):
        try:
            with path.open() as stream:
                rows = list(csv.DictReader(stream))
            values = [int(row["peak_rss_bytes"]) for row in rows]
            if not values or min(values) < 0:
                raise ValueError("invalid RSS samples")
            peaks.append({"file": path.name, "peak_rss_bytes": max(values)})
        except (OSError, ValueError, KeyError, csv.Error) as error:
            errors.append(f"{path.name}: {error}")
    return {"process_peaks": peaks, "errors": errors, "swap_bytes": None,
            "node_peak_bytes": None, "scope": "native per-process high-water RSS; coverage unverified"}


def execute(plan, destination):
    if os.name != "posix":
        raise ValueError("campaign execution currently requires POSIX process groups")
    if any(case["mode"] == "master" for case in plan["manifest"]["cases"]):
        raise ValueError("master campaigns are render-only until queue-worker cancellation is qualified")
    destination = destination.resolve()
    fixture = Path(plan["fixture"])
    if destination == fixture or fixture in destination.parents or destination in fixture.parents:
        raise ValueError("campaign output and fixture must be separate directory trees")
    destination.mkdir(parents=True, exist_ok=True)
    with campaign_lock(destination):
        saved = destination / "campaign.json"
        if saved.exists():
            if json.loads(saved.read_text()) != plan:
                raise ValueError("resume requires identical manifest, fixture and executable fingerprints")
        else:
            with saved.open("x") as stream:
                json.dump(plan, stream, indent=2)
        records = destination / "observations.jsonl"
        history = [json.loads(line) for line in records.read_text().splitlines()] if records.exists() else []
        complete = {r["experiment_id"] for r in history if r["event"] == "finished" and r["status"] == "passed"}
        starts = [r for r in history if r["event"] == "started"]
        budget = plan["manifest"]["budget"]
        failed = False
        for run in plan["runs"]:
            identity = run["experiment_id"]
            if identity in complete:
                continue
            reserved = (len(starts) + 1) * budget["attempt_seconds"]
            if reserved > budget["wall_seconds"] or reserved * budget["cpus"] > budget["cpu_seconds"]:
                raise ValueError("campaign reservation budget exhausted (including prior attempts)")
            if disk_usage(destination) >= budget["disk_bytes"]:
                raise ValueError("campaign disk budget exhausted")
            attempt = sum(r["experiment_id"] == identity for r in starts) + 1
            directory = destination / identity / f"attempt-{attempt}"
            directory.mkdir(parents=True, exist_ok=False)
            start = {"event": "started", "experiment_id": identity, "attempt": attempt,
                     "reserved_seconds": budget["attempt_seconds"], "warmup": run["warmup"]}
            append_record(records, start)
            starts.append(start)
            deadline = time.monotonic() + budget["attempt_seconds"]
            observation = {**start, "event": "finished", "status": "failed",
                           "qualification": "unqualified", "stages": []}
            case = plan["manifest"]["cases"][run["case_index"]]
            try:
                shutil.copytree(fixture, directory / "fixture")
                if tree_files(directory / "fixture") != plan["fixture_files"]:
                    raise ValueError("fixture changed after campaign planning")
                def expand(command):
                    expanded = [s.replace("{workdir}", str(directory)).replace("{fixture}", str(directory / "fixture"))
                                for s in command]
                    expanded[0] = plan["executables"]["hook:" + command[0]]["path"]
                    return expanded
                env = os.environ.copy()
                env.update(OMP_NUM_THREADS=str(run["candidate"]["nthr"]), OMP_DYNAMIC="FALSE",
                           SIMPLE_SEED=str(plan["manifest"]["seed"]))
                key = "private" if case["mode"] == "worker" else "simple"
                commands = [(f"prepare-{i}", expand(command)) for i, command in enumerate(case["prepare"])]
                commands += [("measure", build_command(case, run["candidate"],
                              plan["executables"][key]["path"], directory, directory / "fixture")),
                             ("validate", expand(case["validate"]))]
                for stage, command in commands:
                    if time.monotonic() >= deadline:
                        raise ValueError("attempt budget exhausted before " + stage)
                    if stage == "measure" and not (directory / case["project"]).is_file():
                        raise ValueError("preparation did not create the declared project")
                    outcome = run_process(command, directory, env, deadline,
                                          directory / f"{stage}.stdout", directory / f"{stage}.stderr",
                                          destination, budget["disk_bytes"], plan["cpu_ids"])
                    observation["stages"].append({"stage": stage, "argv": command, **outcome})
                    if outcome["returncode"] != 0 or outcome["reason"]:
                        raise ValueError(f"{stage} failed: {outcome}")
                validation = json.loads((directory / "validate.stdout").read_text())
                canonical(validation)
                if (not isinstance(validation, dict) or validation.get("passed") is not True
                        or not isinstance(validation.get("validator_version"), str)
                        or not validation["validator_version"]
                        or not isinstance(validation.get("metrics"), dict) or not validation["metrics"]):
                    raise ValueError("validator must report passed=true, validator_version and nonempty metrics")
                observation.update(status="passed", validation=validation)
            except (OSError, ValueError, KeyError) as error:
                observation["error"] = str(error)
                failed = True
            finally:
                observation["telemetry"] = telemetry(directory)
                append_record(records, observation)
        return 2 if failed else 0
