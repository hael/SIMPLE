#!/usr/bin/env python3
"""Fit and validate the transparent SIMPLE memory-estimator models."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent.parent
ROOT = SCRIPTS_DIR.parent
sys.path.insert(0, str(SCRIPTS_DIR))

from memory_estimator import DEFAULT_MODELS, round_conservative
def read_ok_rows(paths: Sequence[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as stream:
            rows.extend(row for row in csv.DictReader(stream) if row["status"] == "ok")
    if not rows:
        raise ValueError(f"no successful rows in: {', '.join(map(str, paths))}")
    return rows


def number(row: dict[str, str], name: str) -> float:
    return float(row[name])


FeatureBuilder = Callable[[dict[str, str]], list[float]]


def motion_features(row: dict[str, str]) -> list[float]:
    original_mp = number(row, "pixels_per_frame") / 1.0e6
    effective_mp = number(row, "effective_pixels_per_frame") / 1.0e6
    extra_frames = max(0.0, number(row, "frames") - 4.0)
    extra_threads = max(0.0, number(row, "nthr") - 1.0)
    return [
        original_mp,
        effective_mp,
        original_mp * number(row, "smpd_angstrom_per_pixel"),
        original_mp * extra_frames,
        effective_mp * extra_frames,
        original_mp * extra_threads,
    ]


def ab2d_features(row: dict[str, str]) -> list[float]:
    box = number(row, "box")
    nsample = number(row, "nsample_effective")
    return [
        box**2 / 1.0e4,
        nsample * box**2 / 1.0e6,
        number(row, "ncls"),
        number(row, "nthr"),
    ]


def ab3d_features(row: dict[str, str]) -> list[float]:
    box = number(row, "box")
    smpd = number(row, "smpd_angstrom_per_pixel")
    mask = number(row, "mskdiam_angstrom")
    return [
        box**3 / 1.0e6,
        number(row, "nsample_effective") * box**2 / 1.0e6,
        max(0.0, number(row, "nparts") - 1.0),
        max(0.0, number(row, "nstates") - 1.0),
        (mask / smpd) ** 3 / 1.0e6,
        number(row, "nthr"),
    ]


def design(rows: list[dict[str, str]], builder: FeatureBuilder, target: str) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.asarray([[1.0, *builder(row)] for row in rows], dtype=np.float64)
    observed = np.asarray([number(row, target) for row in rows], dtype=np.float64)
    return matrix, observed


def metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = predicted - observed
    denominator = np.sum((observed - observed.mean()) ** 2)
    return {
        "rmse_mib": float(np.sqrt(np.mean(residual**2))),
        "mae_mib": float(np.mean(np.abs(residual))),
        "mape_percent": float(100.0 * np.mean(np.abs(residual) / observed)),
        "r_squared": float(1.0 - np.sum(residual**2) / denominator),
    }


def leave_one_out(matrix: np.ndarray, observed: np.ndarray) -> np.ndarray:
    predicted = np.empty_like(observed)
    for index in range(len(observed)):
        selected = np.arange(len(observed)) != index
        coefficients = np.linalg.lstsq(matrix[selected], observed[selected], rcond=None)[0]
        predicted[index] = matrix[index] @ coefficients
    return predicted


def conservative_predictions(
    commander: str,
    rows: list[dict[str, str]],
    raw: np.ndarray,
    safety: dict[str, Any],
) -> np.ndarray:
    allocations: list[int] = []
    for row, estimate in zip(rows, raw, strict=True):
        factor = None
        if commander == "motion_correct":
            extra_frames = max(0.0, number(row, "frames") - 4.0)
            factor = float(safety["factor"]) + min(
                float(safety["extra_frame_factor_cap"]),
                float(safety["extra_frame_factor_per_frame"]) * extra_frames,
            )
        allocations.append(round_conservative(float(estimate), safety, factor))
    return np.asarray(allocations, dtype=np.float64)


def fit_one(
    commander: str,
    rows: list[dict[str, str]],
    builder: FeatureBuilder,
    target: str,
    names: list[str],
    model: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    matrix, observed = design(rows, builder, target)
    coefficients = np.linalg.lstsq(matrix, observed, rcond=None)[0]
    raw = matrix @ coefficients
    loo = leave_one_out(matrix, observed)
    allocation = conservative_predictions(commander, rows, raw, model["safety"])
    fitted = deepcopy(model)
    fitted["rows"] = len(rows)
    fitted["coefficients"] = {
        name: float(value) for name, value in zip(["intercept", *names], coefficients, strict=True)
    }
    validation = {
        "commander": commander,
        "rows": len(rows),
        "target": target,
        "fit": metrics(observed, raw),
        "leave_one_out": metrics(observed, loo),
        "allocation_coverage_percent": float(100.0 * np.mean(allocation >= observed)),
        "median_allocation_headroom_mib": float(np.median(allocation - observed)),
        "maximum_underallocation_mib": float(max(0.0, np.max(observed - allocation))),
        "observed_range_mib": [float(observed.min()), float(observed.max())],
        "recommended_range_mib": [int(allocation.min()), int(allocation.max())],
    }
    return fitted, validation


def report_text(validations: list[dict[str, Any]]) -> str:
    lines = [
        "SIMPLE MEMORY ESTIMATOR VALIDATION",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}",
        "",
        "The raw model is ordinary least squares. The recommended allocation adds",
        "a commander-specific safety margin and rounds upward. Coverage below is",
        "in-sample allocation coverage; leave-one-out metrics show model stability.",
        "",
    ]
    for item in validations:
        fit, loo = item["fit"], item["leave_one_out"]
        lines.extend([
            item["commander"],
            f"  Rows: {item['rows']}",
            f"  Target column: {item['target']}",
            f"  Observed range: {item['observed_range_mib'][0]:.1f}-{item['observed_range_mib'][1]:.1f} MiB",
            f"  Raw fit: RMSE {fit['rmse_mib']:.1f} MiB; MAE {fit['mae_mib']:.1f} MiB; R^2 {fit['r_squared']:.3f}",
            f"  Leave-one-out: RMSE {loo['rmse_mib']:.1f} MiB; MAE {loo['mae_mib']:.1f} MiB; R^2 {loo['r_squared']:.3f}",
            f"  Conservative allocation coverage: {item['allocation_coverage_percent']:.1f}%",
            f"  Median allocation headroom: {item['median_allocation_headroom_mib']:.1f} MiB",
            f"  Maximum underallocation: {item['maximum_underallocation_mib']:.1f} MiB",
            f"  Recommended range on calibration data: {item['recommended_range_mib'][0]}-{item['recommended_range_mib'][1]} MiB",
            "",
        ])
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--motion-csv", nargs="+", type=Path, required=True,
        help="one or more motion_correct benchmark results.csv files",
    )
    parser.add_argument("--abinitio2d-csv", type=Path, required=True)
    parser.add_argument("--abinitio3d-csv", type=Path, required=True)
    parser.add_argument("--template", type=Path, default=DEFAULT_MODELS)
    parser.add_argument(
        "--output-models", type=Path, default=ROOT / "output/memory_estimator_models_fitted.json"
    )
    parser.add_argument(
        "--report", type=Path, default=ROOT / "output/memory_estimator_validation.txt"
    )
    parser.add_argument(
        "--install", action="store_true",
        help="replace scripts/memory_estimator_models.json with the fitted model",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    with args.template.open(encoding="utf-8") as stream:
        payload = json.load(stream)
    definitions = [
        ("motion_correct", read_ok_rows(args.motion_csv), motion_features, "peak_rss_mib", [
            "original_megapixels", "effective_megapixels", "original_megapixels_x_smpd",
            "original_megapixels_x_extra_frames", "effective_megapixels_x_extra_frames",
            "original_megapixels_x_extra_threads",
        ]),
        ("abinitio2D", read_ok_rows([args.abinitio2d_csv]), ab2d_features, "peak_rss_mib", [
            "box_squared_per_10000", "sampled_particle_megapixels", "references", "threads",
        ]),
        ("abinitio3D", read_ok_rows([args.abinitio3d_csv]), ab3d_features, "peak_tree_rss_mib", [
            "box_million_voxels", "sampled_particle_megapixels", "extra_partitions",
            "extra_states", "mask_million_voxels", "threads",
        ]),
    ]
    validations: list[dict[str, Any]] = []
    for commander, rows, builder, target, names in definitions:
        fitted, validation = fit_one(
            commander, rows, builder, target, names, payload["models"][commander]
        )
        payload["models"][commander] = fitted
        validations.append(validation)
    payload["generated_utc"] = datetime.now(timezone.utc).isoformat()
    output_models = DEFAULT_MODELS if args.install else args.output_models
    output_models.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    output_models.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.report.write_text(report_text(validations), encoding="utf-8")
    print(output_models)
    print(args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
