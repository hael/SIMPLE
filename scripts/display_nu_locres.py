#!/usr/bin/env python3
"""Display a density map colored by a local-resolution map in ChimeraX.

Recommended: run with ordinary Python. The script checks for the ChimeraX
executable and then relaunches itself inside ChimeraX:

    python3 display_nu_locres.py map.mrc map_nu_locres.mrc

It can also be run directly through ChimeraX:

    chimerax --script "display_nu_locres.py map.mrc map_nu_locres.mrc"

Installation/interface check only:

    python3 display_nu_locres.py --check-install

Options:
    --range BEST_A WORST_A
        Use a fixed color range in Angstrom instead of automatic percentiles.
    --percentiles LOW HIGH
        Percentiles of positive local-resolution voxels used for the automatic
        range (default: 2 98). Ignored when --range is supplied.
    --level VALUE
        Set the density isosurface contour. By default ChimeraX chooses it.
    --palette PALETTE
        ChimeraX palette name or colon-separated colors in increasing-Angstrom
        order (default: blue:cyan:yellow:red).
    --no-key
        Do not draw the color key.
    --no-controls
        Do not open the interactive Surface Color panel. By default this panel
        is opened so color/value changes can be applied to the map with Color.
    --check-install
        Check the ChimeraX executable and required Python interfaces, then exit.

The local-resolution MRC is interpreted as an Angstrom-valued scalar field.
Non-finite and non-positive values are treated as invalid/background when the
automatic range is calculated.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


_IN_CHIMERAX = "session" in globals() and hasattr(globals()["session"], "logger")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Display a density map in ChimeraX and color its isosurface\n"
            "by the Angstrom values in a co-registered local-resolution MRC."
        ),
        epilog=(
            "Examples:\n"
            "  python3 %(prog)s map.mrc map_nu_locres.mrc\n"
            "  python3 %(prog)s map.mrc map_nu_locres.mrc --range 2.6 5.0\n"
            "  python3 %(prog)s map.mrc map_nu_locres.mrc --level 0.08 "
            "--palette ^rdylbu\n"
            "  python3 %(prog)s --check-install"
        ),
    )
    parser.add_argument(
        "density_map",
        nargs="?",
        help="density map whose isosurface will be displayed",
    )
    parser.add_argument(
        "resolution_map",
        nargs="?",
        help="co-registered local-resolution map with voxel values in Angstrom",
    )
    parser.add_argument(
        "--range",
        nargs=2,
        type=float,
        metavar=("BEST_A", "WORST_A"),
        dest="resolution_range",
        help=(
            "fixed color range in Angstrom with BEST_A < WORST_A; values outside "
            "the range use the nearest endpoint color"
        ),
    )
    parser.add_argument(
        "--percentiles",
        nargs=2,
        type=float,
        default=(2.0, 98.0),
        metavar=("LOW", "HIGH"),
        help=(
            "positive-voxel percentiles used to select the automatic color range "
            "(default: 2 98; ignored with --range)"
        ),
    )
    parser.add_argument(
        "--level",
        type=float,
        help="density-map contour level (default: ChimeraX automatic level)",
    )
    parser.add_argument(
        "--palette",
        default="blue:cyan:yellow:red",
        help=(
            "ChimeraX palette name or colon-separated colors in increasing-Angstrom "
            "order (default: blue:cyan:yellow:red; prefix a named palette with ^ "
            "to reverse it)"
        ),
    )
    parser.add_argument(
        "--no-key",
        action="store_true",
        help="do not draw the local-resolution color key",
    )
    parser.add_argument(
        "--no-controls",
        action="store_true",
        help="do not open the interactive Surface Color panel",
    )
    parser.add_argument(
        "--check-install",
        action="store_true",
        help="verify the ChimeraX executable and required Python interfaces, then exit",
    )
    parser.add_argument("--check-interface", action="store_true", help=argparse.SUPPRESS)
    return parser


def _find_chimerax_executable() -> Path | None:
    """Find an executable without importing any ChimeraX Python packages."""
    candidates: list[str] = []
    configured = os.environ.get("CHIMERAX_EXECUTABLE")
    if configured:
        candidates.append(configured)

    for command in ("chimerax", "ChimeraX"):
        found = shutil.which(command)
        if found:
            candidates.append(found)

    # Common locations not necessarily placed on PATH by an installer.
    candidates.extend(
        (
            "/Applications/ChimeraX.app/Contents/MacOS/ChimeraX",
            "/usr/bin/chimerax",
            "/usr/local/bin/chimerax",
        )
    )
    for candidate in candidates:
        path = Path(candidate).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            return path.resolve()
    return None


def _script_argument(argv: list[str]) -> str:
    # ChimeraX --script accepts one string containing the script path and all
    # arguments. shlex quoting preserves spaces and shell metacharacters.
    return shlex.join([str(Path(__file__).resolve()), *argv])


def _launcher_main(argv: list[str]) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    executable = _find_chimerax_executable()
    if executable is None:
        print(
            "ERROR: ChimeraX was not found. Install UCSF ChimeraX and ensure the "
            "'chimerax' executable is on PATH, or set CHIMERAX_EXECUTABLE to its "
            "full path.",
            file=sys.stderr,
        )
        return 2

    if args.check_install:
        command = [
            str(executable),
            "--nogui",
            "--exit",
            "--script",
            _script_argument(["--check-interface"]),
        ]
        print(f"Found ChimeraX executable: {executable}")
    else:
        if args.density_map is None or args.resolution_map is None:
            parser.error(
                "density_map and resolution_map are required unless --check-install is used"
            )
        command = [str(executable), "--script", _script_argument(argv)]

    try:
        return subprocess.call(command)
    except OSError as error:
        print(f"ERROR: Could not start ChimeraX at {executable}: {error}", file=sys.stderr)
        return 2


if not _IN_CHIMERAX:
    raise SystemExit(_launcher_main(sys.argv[1:]))


# These imports must come from ChimeraX's bundled Python environment. Keeping
# them below the launcher branch lets ordinary Python provide a useful install
# error instead of an opaque ModuleNotFoundError.
try:
    import numpy as np

    from chimerax.core.commands import StringArg, run
    from chimerax.core.errors import UserError
except (ImportError, AttributeError) as error:
    print(
        "ERROR: ChimeraX started, but its required Python interface could not be "
        f"loaded ({error}). Repair or update the ChimeraX installation.",
        file=sys.stderr,
    )
    raise SystemExit(2) from error


def _check_chimerax_interface() -> None:
    missing = []
    if not callable(run):
        missing.append("chimerax.core.commands.run")
    if not hasattr(StringArg, "unparse"):
        missing.append("chimerax.core.commands.StringArg.unparse")
    if not hasattr(session, "models"):
        missing.append("session.models")
    if not hasattr(session, "logger"):
        missing.append("session.logger")
    if not hasattr(session, "ui") or not hasattr(session.ui, "is_gui"):
        missing.append("session.ui.is_gui")
    if not hasattr(np, "percentile"):
        missing.append("ChimeraX-bundled numpy.percentile")
    if missing:
        raise UserError("Required ChimeraX interfaces are unavailable: " + ", ".join(missing))


def _require_input_paths(args: argparse.Namespace) -> tuple[str, str]:
    if args.density_map is None or args.resolution_map is None:
        raise UserError("density_map and resolution_map are required; use --help for examples")
    return args.density_map, args.resolution_map


def _existing_map(value: str, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise UserError(f"{label} does not exist or is not a file: {path}")
    recognized = (".mrc", ".mrc.gz", ".map", ".map.gz", ".ccp4", ".ccp4.gz")
    if not path.name.lower().endswith(recognized):
        session.logger.warning(
            f"{label} does not have a usual MRC/CCP4 map suffix: {path.name}. "
            "ChimeraX will still try to identify its format."
        )
    return path


def _first_volume(models, label: str):
    pending = list(models)
    while pending:
        model = pending.pop(0)
        data = getattr(model, "data", None)
        if data is not None and hasattr(data, "matrix"):
            required = ("size", "step", "origin")
            unavailable = [name for name in required if not hasattr(data, name)]
            if unavailable or not hasattr(model, "atomspec"):
                detail = ", ".join(unavailable + ([] if hasattr(model, "atomspec") else ["atomspec"]))
                raise UserError(
                    f"The ChimeraX volume interface for {label} lacks required fields: {detail}"
                )
            return model
        if hasattr(model, "child_models"):
            pending.extend(model.child_models())
    raise UserError(f"Opening {label} did not create a readable ChimeraX volume model")


def _validate_grid_metadata(volume, label: str) -> None:
    try:
        size = tuple(int(value) for value in volume.data.size)
        step = tuple(float(value) for value in volume.data.step)
        origin = tuple(float(value) for value in volume.data.origin)
    except (TypeError, ValueError) as error:
        raise UserError(f"{label} has unreadable grid metadata: {error}") from error
    if len(size) != 3 or any(value < 2 for value in size):
        raise UserError(f"{label} has invalid 3D grid dimensions: {size}")
    if len(step) != 3 or not all(np.isfinite(value) and value > 0.0 for value in step):
        raise UserError(f"{label} has invalid voxel spacing: {step}")
    if len(origin) != 3 or not all(np.isfinite(value) for value in origin):
        raise UserError(f"{label} has invalid grid origin: {origin}")


def _open_volume(path: Path, label: str):
    try:
        models = run(session, f"open {StringArg.unparse(str(path))}")
    except Exception as error:
        raise UserError(f"ChimeraX could not open {label} ({path}): {error}") from error
    volume = _first_volume(models, label)
    _validate_grid_metadata(volume, label)
    return volume


def _matrix_sample(volume, label: str, target_voxels: int = 5_000_000):
    try:
        values = np.asarray(volume.data.matrix())
    except Exception as error:
        raise UserError(f"Could not read voxel values from {label}: {error}") from error
    if values.ndim != 3 or any(length < 2 for length in values.shape):
        raise UserError(f"{label} must be a non-empty 3D grid; observed shape {values.shape}")

    # Limit percentile/statistics work for very large boxes while sampling all
    # three dimensions uniformly. The full map remains loaded and displayed.
    stride = max(1, int(np.ceil((values.size / target_voxels) ** (1.0 / 3.0))))
    sample = values[::stride, ::stride, ::stride]
    return sample, stride


def _inspect_density(volume) -> dict:
    sample, _ = _matrix_sample(volume, "density map")
    finite = sample[np.isfinite(sample)]
    if finite.size == 0:
        raise UserError("density map contains no finite voxel values")
    minimum = float(finite.min())
    maximum = float(finite.max())
    if not minimum < maximum:
        raise UserError("density map is constant and cannot produce a meaningful isosurface")
    nonfinite = int(sample.size - finite.size)
    if nonfinite:
        session.logger.warning(
            f"density map contains {nonfinite:,} non-finite voxels in the sampled grid; "
            "these may cause display artifacts."
        )
    return {"minimum": minimum, "maximum": maximum}


def _inspect_resolution(volume) -> dict:
    sample, stride = _matrix_sample(volume, "resolution map")
    finite_mask = np.isfinite(sample)
    positive = sample[finite_mask & (sample > 0.0)].astype(np.float64, copy=False)
    if positive.size == 0:
        raise UserError(
            "resolution map contains no finite values greater than zero; check that "
            "the supplied volume stores local resolution in Angstrom"
        )

    negative_count = int(np.count_nonzero(finite_mask & (sample < 0.0)))
    nonfinite_count = int(sample.size - np.count_nonzero(finite_mask))
    if negative_count:
        session.logger.warning(
            f"resolution map contains {negative_count:,} negative sampled voxels. They "
            "will be treated as invalid/background for automatic range selection."
        )
    if nonfinite_count:
        session.logger.warning(
            f"resolution map contains {nonfinite_count:,} non-finite sampled voxels. They "
            "will be excluded from automatic range selection."
        )

    p01, p50, p99 = (float(value) for value in np.percentile(positive, (1.0, 50.0, 99.0)))
    if p01 < 0.5 or p99 > 100.0:
        session.logger.warning(
            f"The sampled local-resolution range is unusual ({p01:.3g}-{p99:.3g} "
            "Angstrom at the 1st-99th percentiles). Verify that the resolution map "
            "contains values in Angstrom rather than density or reciprocal resolution."
        )

    voxel_spacing = max(abs(float(value)) for value in volume.data.step)
    nyquist = 2.0 * voxel_spacing
    if p01 < 0.95 * nyquist:
        session.logger.warning(
            f"Some local-resolution values ({p01:.3f} Angstrom at the 1st percentile) "
            f"are better than the approximate map Nyquist limit ({nyquist:.3f} "
            "Angstrom). Check the MRC voxel spacing and map provenance."
        )

    return {
        "positive": positive,
        "minimum": float(positive.min()),
        "maximum": float(positive.max()),
        "median": p50,
        "stride": stride,
    }


def _auto_range(positive, percentiles: tuple[float, float]) -> tuple[float, float]:
    low_p, high_p = percentiles
    if not (0.0 <= low_p < high_p <= 100.0):
        raise UserError("--percentiles must satisfy 0 <= LOW < HIGH <= 100")

    low, high = (float(value) for value in np.percentile(positive, (low_p, high_p)))
    if not low < high:
        low = float(positive.min())
        high = float(positive.max())
    if not low < high:
        # A constant local-resolution map can still be displayed.
        pad = max(abs(low) * 0.01, 0.01)
        low -= pad
        high += pad
    return low, high


def _warn_if_grids_differ(density, resolution) -> None:
    checks = [
        ("grid dimensions", tuple(density.data.size), tuple(resolution.data.size)),
        ("voxel spacing", tuple(density.data.step), tuple(resolution.data.step)),
        ("origin", tuple(density.data.origin), tuple(resolution.data.origin)),
    ]
    if hasattr(density.data, "rotation") and hasattr(resolution.data, "rotation"):
        checks.append(("grid rotation", density.data.rotation, resolution.data.rotation))
    differences = [name for name, left, right in checks if not np.allclose(left, right)]
    if differences:
        session.logger.warning(
            "The density and local-resolution maps differ in "
            + ", ".join(differences)
            + ". ChimeraX will interpolate by physical coordinates, but verify that "
            "the two maps are correctly registered."
        )


def _validate_requested_range(low: float, high: float, stats: dict) -> None:
    if not (np.isfinite(low) and np.isfinite(high) and 0.0 < low < high):
        raise UserError("--range must contain finite positive values with BEST_A < WORST_A")
    if high < stats["minimum"] or low > stats["maximum"]:
        raise UserError(
            f"--range {low:g}-{high:g} Angstrom does not overlap the positive "
            f"resolution-map values ({stats['minimum']:.3g}-{stats['maximum']:.3g} "
            "Angstrom)"
        )


def main(argv: list[str]) -> None:
    _check_chimerax_interface()
    args = _parser().parse_args(argv)

    if args.check_interface or args.check_install:
        run(session, "version")
        session.logger.info(
            "ChimeraX installation check passed: executable, session, command runner, "
            "string parser, model manager, logger, UI interface, and bundled NumPy "
            "are available."
        )
        return

    density_input, resolution_input = _require_input_paths(args)
    density_path = _existing_map(density_input, "density map")
    resolution_path = _existing_map(resolution_input, "resolution map")
    if density_path == resolution_path:
        raise UserError("density_map and resolution_map resolve to the same file")

    density = _open_volume(density_path, "density map")
    resolution = _open_volume(resolution_path, "resolution map")
    density_spec = density.atomspec
    resolution_spec = resolution.atomspec

    _warn_if_grids_differ(density, resolution)
    density_stats = _inspect_density(density)
    resolution_stats = _inspect_resolution(resolution)

    if args.resolution_range is None:
        low, high = _auto_range(resolution_stats["positive"], tuple(args.percentiles))
        range_source = (
            f"automatic {args.percentiles[0]:g}-{args.percentiles[1]:g} percentile "
            f"range from a stride-{resolution_stats['stride']} sample of positive voxels"
        )
    else:
        low, high = args.resolution_range
        _validate_requested_range(low, high, resolution_stats)
        range_source = "user-specified range"

    if args.level is not None:
        if not np.isfinite(args.level):
            raise UserError("--level must be finite")
        if not density_stats["minimum"] < args.level < density_stats["maximum"]:
            raise UserError(
                f"--level {args.level:g} is outside the sampled density range "
                f"({density_stats['minimum']:.6g}-{density_stats['maximum']:.6g})"
            )

    # Show only the density isosurface. The local-resolution map remains open as
    # the scalar sampling source, but its own non-density surface is hidden.
    run(session, f"volume {density_spec} style surface show")
    run(session, f"volume {resolution_spec} hide")
    if args.level is not None:
        run(session, f"volume {density_spec} level {args.level:.9g}")

    palette = StringArg.unparse(args.palette)
    key = "false" if args.no_key else "true"
    run(
        session,
        f"color sample {density_spec} map {resolution_spec} "
        f"palette {palette} range {low:.6g},{high:.6g} key {key}",
    )
    run(session, f"view {density_spec}")

    if session.ui.is_gui:
        # `key true` creates the correct legend but also raises the Color Key
        # editor. That editor controls only the legend and cannot recolor the
        # map, which is easy to misinterpret. Keep the legend visible, hide its
        # editor, and show the Surface Color panel that actually reapplies
        # palette/value changes to the density surface.
        if not args.no_key:
            run(session, "ui tool hide 'Color Key'", log=False)
        if not args.no_controls:
            run(session, "ui tool show 'Surface Color'", log=False)
            session.logger.info(
                "Interactive recoloring: edit colors or Angstrom levels in the "
                "Surface Color panel, click Color to update the density surface, "
                "then click Key to synchronize the legend. The separate Color Key "
                "tool edits only the legend."
            )

    session.logger.info(
        f"Displayed {density_path.name} colored by {resolution_path.name}. "
        f"Local-resolution color range: {low:.3f}-{high:.3f} Angstrom "
        f"({range_source}); sampled median {resolution_stats['median']:.3f} Angstrom. "
        "Lower values are better resolution. Non-positive values were treated as "
        "background for automatic range selection."
    )


main(sys.argv[1:])
