#!/usr/bin/env python3
"""Summarize cluster_cavgs_quality quality_mode=analyze outputs."""

from __future__ import annotations

import argparse
import csv
import io
import math
from pathlib import Path
from typing import Any


FEATURES = "cavgs_quality_features.txt"
SUMMARY = "cavgs_quality_reference_summary.txt"
THRESHOLDS = "cavgs_quality_reference_threshold_scan.txt"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize SIMPLE cluster_cavgs_quality analyze-mode diagnostics."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        default=["."],
        help="Run directories, output files, or parent directories when --recursive is set.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Find run directories below each input path by searching for cavgs_quality_features.txt.",
    )
    parser.add_argument(
        "--top-classes",
        type=int,
        default=8,
        help="Maximum false-positive and false-negative classes to print per run.",
    )
    parser.add_argument(
        "--top-features",
        type=int,
        default=8,
        help="Maximum features to print in ranked feature summaries.",
    )
    return parser.parse_args()


def as_float(value: Any, default: float = math.nan) -> float:
    text = str(value).strip()
    if not text:
        return default
    try:
        return float(text)
    except ValueError:
        return default


def as_int(value: Any, default: int = 0) -> int:
    text = str(value).strip()
    if not text:
        return default
    try:
        return int(float(text))
    except ValueError:
        return default


def as_bool(value: Any) -> bool:
    return str(value).strip().upper() in {"T", "TRUE", "1", "YES", ".TRUE."}


def fmt(value: Any, digits: int = 3) -> str:
    number = as_float(value)
    if math.isnan(number):
        return "n/a"
    if abs(number) >= 1000 or (0 < abs(number) < 0.001):
        return f"{number:.{digits}e}"
    return f"{number:.{digits}f}"


def clean_row(row: dict[str, Any]) -> dict[str, str]:
    return {str(key).strip(): str(value).strip() for key, value in row.items() if key is not None}


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return [clean_row(row) for row in csv.DictReader(handle)]


def read_summary(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    if not path.exists():
        return {}, []
    lines = path.read_text().splitlines()
    table_start = None
    values: dict[str, str] = {}
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("feature,"):
            table_start = index
            break
        if not stripped or stripped.startswith("#"):
            continue
        if "=" in stripped:
            key, value = stripped.split("=", 1)
            values[key.strip()] = value.strip()
    feature_rows: list[dict[str, str]] = []
    if table_start is not None:
        text = "\n".join(lines[table_start:])
        feature_rows = [clean_row(row) for row in csv.DictReader(io.StringIO(text))]
    return values, feature_rows


def expand_run_dirs(paths: list[str], recursive: bool) -> list[Path]:
    run_dirs: set[Path] = set()
    for raw in paths:
        path = Path(raw).expanduser()
        if path.is_file():
            run_dirs.add(path.parent.resolve())
        elif recursive:
            for found in path.rglob(FEATURES):
                run_dirs.add(found.parent.resolve())
        else:
            run_dirs.add(path.resolve())
    return sorted(run_dirs)


def confusion_from_rows(rows: list[dict[str, str]]) -> dict[str, int]:
    counts = {"tp": 0, "fp": 0, "tn": 0, "fn": 0}
    for row in rows:
        if "manual_state" not in row:
            continue
        auto = as_int(row.get("state")) > 0
        manual = as_int(row.get("manual_state")) > 0
        if auto and manual:
            counts["tp"] += 1
        elif auto and not manual:
            counts["fp"] += 1
        elif not auto and not manual:
            counts["tn"] += 1
        else:
            counts["fn"] += 1
    return counts


def binary_metrics(tp: int, fp: int, tn: int, fn: int) -> dict[str, float]:
    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    specificity = safe_div(tn, tn + fp)
    f1 = safe_div(2 * precision * recall, precision + recall)
    balacc = 0.5 * (recall + specificity)
    accuracy = safe_div(tp + tn, tp + fp + tn + fn)
    return {
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "f1": f1,
        "balanced_accuracy": balacc,
        "accuracy": accuracy,
    }


def safe_div(num: float, den: float) -> float:
    return num / den if den else math.nan


def weights_from_features(feature_rows: list[dict[str, str]]) -> dict[str, float]:
    weights: dict[str, float] = {}
    for row in feature_rows:
        name = row.get("feature", "")
        if name:
            weights[name] = as_float(row.get("current_weight"), 0.0)
    return weights


def z_feature_names(rows: list[dict[str, str]]) -> list[str]:
    if not rows:
        return []
    return sorted(key[2:] for key in rows[0] if key.startswith("z_"))


def feature_signal(
    row: dict[str, str],
    features: list[str],
    weights: dict[str, float],
    strongest: bool,
    limit: int = 3,
) -> str:
    scored: list[tuple[float, str, float, float]] = []
    use_weights = any(abs(weight) > 0 for weight in weights.values())
    for feature in features:
        zval = as_float(row.get(f"z_{feature}"))
        if math.isnan(zval):
            continue
        weight = weights.get(feature, 1.0 if not use_weights else 0.0)
        score = zval * weight if use_weights else zval
        scored.append((score, feature, zval, weight))
    scored.sort(reverse=strongest, key=lambda item: item[0])
    parts = []
    for _, feature, zval, weight in scored[:limit]:
        if use_weights:
            parts.append(f"{feature}:z={fmt(zval)} w={fmt(weight, 2)}")
        else:
            parts.append(f"{feature}:z={fmt(zval)}")
    return "; ".join(parts) if parts else "n/a"


def class_line(row: dict[str, str], features: list[str], weights: dict[str, float], strongest: bool) -> str:
    cls = row.get("class", "?")
    score = fmt(row.get("quality_score"))
    cluster = row.get("quality_cluster", "?")
    hard = "T" if as_bool(row.get("hard_reject")) else "F"
    signal = feature_signal(row, features, weights, strongest=strongest)
    return f"- class {cls}: score={score}, cluster={cluster}, hard_reject={hard}; {signal}"


def best_threshold_row(rows: list[dict[str, str]], metric: str) -> dict[str, str] | None:
    if not rows:
        return None
    return max(rows, key=lambda row: as_float(row.get(metric), -math.inf))


def nearest_threshold_row(rows: list[dict[str, str]], threshold: float) -> dict[str, str] | None:
    if not rows or math.isnan(threshold):
        return None
    return min(rows, key=lambda row: abs(as_float(row.get("score_threshold")) - threshold))


def summarize_feature_rows(feature_rows: list[dict[str, str]], limit: int) -> list[str]:
    if not feature_rows:
        return ["- Feature summary missing."]
    ranked = sorted(
        feature_rows,
        key=lambda row: (as_float(row.get("auc"), -math.inf), abs(as_float(row.get("robust_separation"), 0.0))),
        reverse=True,
    )
    lines = ["- Strongest features by AUC/separation:"]
    for row in ranked[:limit]:
        lines.append(
            "  "
            + f"{row.get('feature', '?')}: auc={fmt(row.get('auc'))}, "
            + f"sep={fmt(row.get('robust_separation'))}, "
            + f"weight={fmt(row.get('current_weight'), 2)}, "
            + f"suggested={fmt(row.get('suggested_weight'), 2)}"
        )
    weak = [
        row
        for row in feature_rows
        if as_float(row.get("auc"), 1.0) < 0.5 and as_float(row.get("current_weight"), 0.0) > 0
    ]
    if weak:
        weak.sort(key=lambda row: as_float(row.get("auc"), 1.0))
        lines.append("- Weighted features with AUC below 0.5:")
        for row in weak[:limit]:
            lines.append(
                "  "
                + f"{row.get('feature', '?')}: auc={fmt(row.get('auc'))}, "
                + f"sep={fmt(row.get('robust_separation'))}, weight={fmt(row.get('current_weight'), 2)}"
            )
    return lines


def summarize_thresholds(
    threshold_rows: list[dict[str, str]], summary: dict[str, str]
) -> list[str]:
    if not threshold_rows:
        return ["- Threshold scan missing."]
    lines: list[str] = []
    current_threshold = as_float(summary.get("current_score_threshold"))
    current = nearest_threshold_row(threshold_rows, current_threshold)
    best_bal = best_threshold_row(threshold_rows, "balanced_accuracy")
    best_f1 = best_threshold_row(threshold_rows, "f1")
    if current:
        lines.append(
            "- Current threshold: "
            + f"{fmt(current.get('score_threshold'))}, selected={current.get('selected', '?')}, "
            + f"F1={fmt(current.get('f1'))}, balacc={fmt(current.get('balanced_accuracy'))}"
        )
    if best_bal:
        lines.append(
            "- Best balanced-accuracy threshold: "
            + f"{fmt(best_bal.get('score_threshold'))}, selected={best_bal.get('selected', '?')}, "
            + f"F1={fmt(best_bal.get('f1'))}, balacc={fmt(best_bal.get('balanced_accuracy'))}"
        )
    if best_f1 and best_f1 is not best_bal:
        lines.append(
            "- Best F1 threshold: "
            + f"{fmt(best_f1.get('score_threshold'))}, selected={best_f1.get('selected', '?')}, "
            + f"F1={fmt(best_f1.get('f1'))}, balacc={fmt(best_f1.get('balanced_accuracy'))}"
        )
    return lines


def summarize_run(run_dir: Path, top_classes: int, top_features: int) -> tuple[str, dict[str, Any]]:
    feature_path = run_dir / FEATURES
    summary_path = run_dir / SUMMARY
    threshold_path = run_dir / THRESHOLDS

    rows = read_csv(feature_path)
    summary, feature_rows = read_summary(summary_path)
    threshold_rows = read_csv(threshold_path)
    weights = weights_from_features(feature_rows)
    z_names = z_feature_names(rows)

    lines = [f"## {run_dir}"]
    if not rows:
        lines.append(f"- Missing `{FEATURES}`.")
        return "\n".join(lines), {"run": str(run_dir), "ok": False}

    counts = {
        "tp": as_int(summary.get("true_positive")),
        "fp": as_int(summary.get("false_positive")),
        "tn": as_int(summary.get("true_negative")),
        "fn": as_int(summary.get("false_negative")),
    }
    if sum(counts.values()) == 0:
        counts = confusion_from_rows(rows)
    metrics = binary_metrics(counts["tp"], counts["fp"], counts["tn"], counts["fn"])
    n_classes = as_int(summary.get("n_classes"), len(rows))
    manual_selected = as_int(
        summary.get("manual_selected"),
        sum(1 for row in rows if as_int(row.get("manual_state")) > 0),
    )
    auto_selected = as_int(
        summary.get("auto_selected"),
        sum(1 for row in rows if as_int(row.get("state")) > 0),
    )
    hard_manual_good = sum(
        1
        for row in rows
        if as_int(row.get("manual_state")) > 0 and as_bool(row.get("hard_reject"))
    )

    lines.append(
        "- Overall: "
        + f"n={n_classes}, manual_selected={manual_selected}, auto_selected={auto_selected}, "
        + f"TP/FP/TN/FN={counts['tp']}/{counts['fp']}/{counts['tn']}/{counts['fn']}"
    )
    lines.append(
        "- Metrics: "
        + f"precision={fmt(summary.get('precision', metrics['precision']))}, "
        + f"recall={fmt(summary.get('recall', metrics['recall']))}, "
        + f"specificity={fmt(summary.get('specificity', metrics['specificity']))}, "
        + f"F1={fmt(summary.get('f1', metrics['f1']))}, "
        + f"balacc={fmt(summary.get('balanced_accuracy', metrics['balanced_accuracy']))}, "
        + f"AUC={fmt(summary.get('score_auc'))}"
    )
    lines.append(
        "- Thresholds: "
        + f"raw={fmt(summary.get('raw_score_threshold'))}, "
        + f"margin={fmt(summary.get('threshold_boundary_margin'))}, "
        + f"effective={fmt(summary.get('current_score_threshold'))}, "
        + f"hard_rejected_manual_good={hard_manual_good}"
    )
    lines.extend(summarize_thresholds(threshold_rows, summary))
    lines.extend(summarize_feature_rows(feature_rows, top_features))

    false_negatives = [
        row
        for row in rows
        if as_int(row.get("manual_state")) > 0 and as_int(row.get("state")) <= 0
    ]
    false_positives = [
        row
        for row in rows
        if as_int(row.get("manual_state")) <= 0 and as_int(row.get("state")) > 0
    ]
    false_negatives.sort(key=lambda row: as_float(row.get("quality_score")), reverse=True)
    false_positives.sort(key=lambda row: as_float(row.get("quality_score")), reverse=True)

    if false_negatives:
        lines.append(f"- False negatives, highest scores first (showing {min(top_classes, len(false_negatives))}):")
        for row in false_negatives[:top_classes]:
            lines.append(class_line(row, z_names, weights, strongest=False))
    if false_positives:
        lines.append(f"- False positives, highest scores first (showing {min(top_classes, len(false_positives))}):")
        for row in false_positives[:top_classes]:
            lines.append(class_line(row, z_names, weights, strongest=True))
    if not false_negatives and not false_positives:
        lines.append("- Auto selection matches the manual reference for all classes in the feature table.")

    record = {
        "run": str(run_dir),
        "ok": True,
        "n": n_classes,
        "manual_selected": manual_selected,
        "auto_selected": auto_selected,
        **counts,
        "precision": as_float(summary.get("precision"), metrics["precision"]),
        "recall": as_float(summary.get("recall"), metrics["recall"]),
        "specificity": as_float(summary.get("specificity"), metrics["specificity"]),
        "f1": as_float(summary.get("f1"), metrics["f1"]),
        "balanced_accuracy": as_float(summary.get("balanced_accuracy"), metrics["balanced_accuracy"]),
        "auc": as_float(summary.get("score_auc")),
        "hard_manual_good": hard_manual_good,
    }
    return "\n".join(lines), record


def comparison_table(records: list[dict[str, Any]]) -> str:
    ok_records = [record for record in records if record.get("ok")]
    if len(ok_records) < 2:
        return ""
    lines = [
        "## Cross-run comparison",
        "| run | n | manual | auto | TP/FP/TN/FN | precision | recall | F1 | balacc | AUC | hard manual-good |",
        "| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for record in ok_records:
        run_name = Path(str(record["run"])).name or str(record["run"])
        lines.append(
            f"| {run_name} | {record['n']} | {record['manual_selected']} | {record['auto_selected']} | "
            + f"{record['tp']}/{record['fp']}/{record['tn']}/{record['fn']} | "
            + f"{fmt(record['precision'])} | {fmt(record['recall'])} | {fmt(record['f1'])} | "
            + f"{fmt(record['balanced_accuracy'])} | {fmt(record['auc'])} | {record['hard_manual_good']} |"
        )
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    run_dirs = expand_run_dirs(args.paths, args.recursive)
    if not run_dirs:
        print("No run directories found.")
        return 1
    records: list[dict[str, Any]] = []
    blocks: list[str] = []
    for run_dir in run_dirs:
        block, record = summarize_run(run_dir, args.top_classes, args.top_features)
        blocks.append(block)
        records.append(record)
    table = comparison_table(records)
    if table:
        blocks.append(table)
    print("\n\n".join(blocks))
    return 0 if any(record.get("ok") for record in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
