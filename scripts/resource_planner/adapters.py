"""Commander capabilities and argv construction; no process launching here."""

import re


CAPABILITIES = {
    "pick": {"modes": ("shared",), "unit": "micrograph"},
    "solve3D": {"modes": ("shared",), "unit": "particle"},
    "solve2D": {"modes": ("shared",), "unit": "particle"},
    "motion_correct": {"modes": ("worker", "master"), "unit": "movie"},
    "ctf_estimate": {"modes": ("worker", "master"), "unit": "micrograph"},
    "ctf_correct": {"modes": ("shared",), "unit": "particle"},
}
OWNED = {"prg", "projfile", "mkdir", "nthr", "nparts", "ncunits", "part",
         "fromp", "top", "qsys_name", "memreport", "memreport_interval",
         "outfile", "outstk", "stream", "job_memory_per_task"}


def positive_int(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def validate_case(case, cpus):
    name = case["commander"]
    if name not in CAPABILITIES or case["mode"] not in CAPABILITIES[name]["modes"]:
        raise ValueError(f"unsupported commander/mode: {name}/{case['mode']}")
    positive_int(case["work_objects"], "work_objects")
    science = case["science"]
    if not isinstance(science, dict) or not science:
        raise ValueError("science must explicitly specify the scientific settings")
    for key, value in science.items():
        if not re.fullmatch(r"[a-z][a-z0-9_]*", key) or key in OWNED:
            raise ValueError(f"invalid or runner-owned science key: {key}")
        if type(value) not in (str, int, float) or not str(value) or '\x00' in str(value):
            raise ValueError(f"invalid value for {key}")
    required = {"pick": ("picker", "moldiam", "pcontrast", "lp"),
                "solve3D": ("pgrp", "nstates", "nstages", "nsample", "mskdiam", "vol1"),
                "solve2D": ("ncls", "mskdiam", "nstages", "nits_per_stage", "extr_lim",
                              "eo_stage", "autoscale", "refine", "nsample", "ctf"),
                "motion_correct": ("algorithm", "mcpatch", "downscale"),
                "ctf_estimate": ("ctfpatch", "pspecsz", "fit_phshift", "dfmin", "dfmax", "hp", "lp"),
                "ctf_correct": ("ctf_correct_mode",)}[name]
    if any(key not in science for key in required):
        raise ValueError(f"{name} requires explicit science keys: {', '.join(required)}")
    if name == "ctf_correct" and science["ctf_correct_mode"] not in ("phaseflip", "wiener"):
        raise ValueError("ctf_correct_mode must be phaseflip or wiener")
    candidates = case["candidates"]
    if not isinstance(candidates, list) or not candidates:
        raise ValueError("candidates must be a nonempty array")
    for candidate in candidates:
        expected = {"nthr"}
        if case["mode"] == "master":
            expected |= {"nparts", "ncunits", "job_memory_per_task"}
        if set(candidate) != expected:
            raise ValueError(f"{name}/{case['mode']} candidate requires exactly {sorted(expected)}")
        for key, value in candidate.items():
            positive_int(value, key)
        if candidate["nthr"] * candidate.get("ncunits", 1) > cpus:
            raise ValueError("candidate exceeds allocated CPU budget")
        if case["mode"] == "master" and not candidate["ncunits"] <= candidate["nparts"] <= case["work_objects"]:
            raise ValueError("requires ncunits <= nparts <= work_objects")


def build_command(case, candidate, executable, workdir, fixture):
    """Master argv intentionally omits worker selectors; worker argv supplies all."""
    def expand(value):
        return str(value).replace("{workdir}", str(workdir)).replace("{fixture}", str(fixture))
    params = {**case["science"], **candidate, "prg": case["commander"],
              "projfile": str(workdir / case["project"]), "mkdir": "no",
              "qsys_name": "local", "memreport": "yes", "memreport_interval": 1}
    if case["mode"] == "worker":
        params.update(part=1, nparts=1, fromp=1, top=case["work_objects"],
                      outfile="benchmark_oris.simple")
    elif case["commander"] == "pick":
        params.update(fromp=1, top=case["work_objects"], nparts=1, part=1, outfile="benchmark_oris.simple")
    elif case["mode"] == "shared" and case["commander"] == "ctf_correct":
        params["outstk"] = "benchmark_corrected.mrcs"
    # SIMPLE's parser requires the program selector before every other argument.
    program = params.pop("prg")
    return [str(executable), f"prg={program}"] + [
        f"{key}={expand(value)}" for key, value in sorted(params.items())]
