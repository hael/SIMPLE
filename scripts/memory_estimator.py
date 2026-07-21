#!/usr/bin/env python3
"""Estimate peak memory for calibrated SIMPLE commanders.

The regression estimates measured peak RSS.  The reported recommendation adds
an allocation safety margin and rounds upward to a scheduler-friendly block.
Inputs outside the calibration domain are accepted, but explicitly flagged as
extrapolations.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable, Sequence


DEFAULT_MODELS = Path(__file__).with_name("memory_estimator_models.json")


def positive_int(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return value


def nonnegative_int(text: str) -> int:
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return value


def positive_float(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or value <= 0.0:
        raise argparse.ArgumentTypeError("must be finite and greater than zero")
    return value


def load_models(path: Path = DEFAULT_MODELS) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        payload = json.load(stream)
    if payload.get("schema_version") != 1:
        raise ValueError(f"unsupported model schema in {path}")
    return payload


def round_conservative(raw_mib: float, safety: dict[str, Any], factor: float | None = None) -> int:
    multiplier = float(safety["factor"] if factor is None else factor)
    protected = max(raw_mib * multiplier, raw_mib + float(safety["floor_mib"]))
    block = int(safety["round_mib"])
    return max(block, math.ceil(protected / block) * block)


def effective_dimension(size: int, smpd: float, target_smpd: float) -> int:
    scale = min(1.0, smpd / target_smpd)
    if scale >= 0.99:
        return size
    value = size * scale
    rounded = math.floor(value + 0.5)
    if rounded % 2 == 0:
        return rounded
    return rounded - 1 if abs(rounded - 1 - value) <= abs(rounded + 1 - value) else rounded + 1


def range_warnings(inputs: dict[str, float], calibration: dict[str, Any]) -> list[str]:
    warnings: list[str] = []
    for name, value in inputs.items():
        bounds = calibration.get(name)
        if isinstance(bounds, list) and len(bounds) == 2 and not bounds[0] <= value <= bounds[1]:
            warnings.append(
                f"{name}={value:g} is outside calibrated range [{bounds[0]:g}, {bounds[1]:g}]"
            )
    return warnings


def estimate_motion_correct(args: argparse.Namespace, model: dict[str, Any]) -> dict[str, Any]:
    c = model["coefficients"]
    original_mp = args.xdim * args.ydim / 1.0e6
    effective_x = effective_dimension(args.xdim, args.smpd, args.smpd_downscale)
    effective_y = effective_dimension(args.ydim, args.smpd, args.smpd_downscale)
    effective_mp = effective_x * effective_y / 1.0e6
    extra_frames = max(0, args.frames - 4)
    extra_threads = max(0, args.threads - 1)
    raw = (
        c["intercept"]
        + c["original_megapixels"] * original_mp
        + c["effective_megapixels"] * effective_mp
        + c["original_megapixels_x_smpd"] * original_mp * args.smpd
        + c["original_megapixels_x_extra_frames"] * original_mp * extra_frames
        + c["effective_megapixels_x_extra_frames"] * effective_mp * extra_frames
        + c["original_megapixels_x_extra_threads"] * original_mp * extra_threads
    )
    safety = model["safety"]
    factor = float(safety["factor"]) + min(
        float(safety["extra_frame_factor_cap"]),
        float(safety["extra_frame_factor_per_frame"]) * extra_frames,
    )
    warnings = range_warnings(
        {"xdim": args.xdim, "ydim": args.ydim, "smpd": args.smpd,
         "frames": args.frames, "threads": args.threads,
         "smpd_downscale": args.smpd_downscale},
        model["calibration"],
    )
    if args.xdim != args.ydim:
        warnings.append("rectangular movies were not present in the calibration data")
    return result_payload("motion_correct", model, raw, factor, warnings, {
        "effective_xdim": effective_x, "effective_ydim": effective_y,
    })


def sampled_particles(requested: int, total: int) -> int:
    return total if requested == 0 else min(requested, total)


def estimate_abinitio2d(args: argparse.Namespace, model: dict[str, Any]) -> dict[str, Any]:
    c = model["coefficients"]
    nsample = sampled_particles(args.sampled_particles, args.particles)
    raw = (
        c["intercept"]
        + c["box_squared_per_10000"] * args.box**2 / 1.0e4
        + c["sampled_particle_megapixels"] * nsample * args.box**2 / 1.0e6
        + c["references"] * args.references
        + c["threads"] * args.threads
    )
    warnings = range_warnings(
        {"particles": args.particles, "box": args.box, "threads": args.threads,
         "references": args.references}, model["calibration"]
    )
    return result_payload("abinitio2D", model, raw, None, warnings, {
        "sampled_particles_effective": nsample,
    })


def estimate_abinitio3d(args: argparse.Namespace, model: dict[str, Any]) -> dict[str, Any]:
    c = model["coefficients"]
    nsample = sampled_particles(args.sampled_particles, args.particles)
    mask_mvox = (args.mask_diameter / args.smpd) ** 3 / 1.0e6
    raw = (
        c["intercept"]
        + c["box_million_voxels"] * args.box**3 / 1.0e6
        + c["sampled_particle_megapixels"] * nsample * args.box**2 / 1.0e6
        + c["extra_partitions"] * max(0, args.partitions - 1)
        + c["extra_states"] * max(0, args.states - 1)
        + c["mask_million_voxels"] * mask_mvox
        + c["threads"] * args.threads
        + c.get("capped_extra_iterations", 0.0) * 4
    )
    warnings = range_warnings(
        {"particles": args.particles, "box": args.box, "smpd": args.smpd,
         "mask_diameter": args.mask_diameter, "threads": args.threads,
         "partitions": args.partitions, "states": args.states,
         "iterations": 20}, model["calibration"]
    )
    return result_payload("abinitio3D", model, raw, None, warnings, {
        "sampled_particles_effective": nsample,
        "stage1_max_iterations": 20,
    })


def result_payload(
    commander: str,
    model: dict[str, Any],
    raw_mib: float,
    factor: float | None,
    warnings: list[str],
    derived: dict[str, Any],
) -> dict[str, Any]:
    return {
        "commander": commander,
        "target": model["target"],
        "raw_estimate_mib": round(raw_mib, 3),
        "recommended_memory_mib": round_conservative(raw_mib, model["safety"], factor),
        "calibrated": not warnings,
        "warnings": warnings,
        "derived": derived,
    }


def add_common_sample_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--sampled-particles", type=nonnegative_int, default=0,
        help="0 means all particles",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    sub = parser.add_subparsers(dest="commander", required=True)

    motion = sub.add_parser("motion_correct")
    motion.add_argument("--xdim", type=positive_int, required=True)
    motion.add_argument("--ydim", type=positive_int, required=True)
    motion.add_argument("--frames", type=positive_int, required=True)
    motion.add_argument("--smpd", type=positive_float, required=True)
    motion.add_argument("--smpd-downscale", type=positive_float, default=1.3)
    motion.add_argument("--threads", type=positive_int, default=1)

    ab2d = sub.add_parser("abinitio2D")
    ab2d.add_argument("--particles", type=positive_int, required=True)
    ab2d.add_argument("--box", type=positive_int, required=True)
    ab2d.add_argument("--threads", type=positive_int, default=4)
    ab2d.add_argument("--references", type=positive_int, required=True)
    add_common_sample_argument(ab2d)

    ab3d = sub.add_parser("abinitio3D")
    ab3d.add_argument("--particles", type=positive_int, required=True)
    ab3d.add_argument("--box", type=positive_int, required=True)
    ab3d.add_argument("--smpd", type=positive_float, required=True)
    ab3d.add_argument("--mask-diameter", type=positive_float, required=True)
    ab3d.add_argument("--threads", type=positive_int, default=4)
    ab3d.add_argument("--partitions", type=positive_int, default=1)
    ab3d.add_argument("--states", type=positive_int, default=1)
    add_common_sample_argument(ab3d)
    return parser


ESTIMATORS: dict[str, Callable[[argparse.Namespace, dict[str, Any]], dict[str, Any]]] = {
    "motion_correct": estimate_motion_correct,
    "abinitio2D": estimate_abinitio2d,
    "abinitio3D": estimate_abinitio3d,
}


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    payload = load_models(args.models)
    model = payload["models"][args.commander]
    result = ESTIMATORS[args.commander](args, model)
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"{result['commander']}: {result['recommended_memory_mib']} MiB recommended")
        print(f"Raw fitted peak: {result['raw_estimate_mib']:.1f} MiB ({result['target']})")
        for warning in result["warnings"]:
            print(f"WARNING: {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
