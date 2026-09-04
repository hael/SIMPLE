"""Read-only 2D class-selection data extracted from a SIMPLE project."""

import hashlib
import math
import os
import struct
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .mrc import read_mrc_stack_info


class ClassSelectionError(ValueError):
    """Raised when a project cannot supply a safe 2D class-selection dataset."""


@dataclass(frozen=True)
class SIMPLEProjectSegment:
    """One fixed-width segment described by the SIMPLE project header."""

    index: int
    oritype: str
    fromp: int
    top: int
    record_width: int
    record_count: int
    first_data_byte: int


class SIMPLEProjectFileReader:
    """Read the small subset of a binary SIMPLE project needed by the viewer."""

    MAX_SEGMENTS = 20
    HEADER_VALUES_PER_SEGMENT = 5
    INT64_SIZE = 8
    HEADER_SIZE = MAX_SEGMENTS * HEADER_VALUES_PER_SEGMENT * INT64_SIZE
    PARTICLE_SEGMENTS = frozenset((3, 6))
    SEGMENT_NAMES = {
        1: "mic",
        2: "stk",
        3: "ptcl2D",
        4: "cls2D",
        5: "cls3D",
        6: "ptcl3D",
        7: "out",
        8: "optics",
        9: "projinfo",
        10: "jobproc",
        11: "compenv",
    }

    def __init__(self, project_path):
        self.project_path = str(project_path)
        self._segments = None

    def _header_values(self):
        with open(self.project_path, "rb") as project_file:
            header = project_file.read(self.HEADER_SIZE)
        if len(header) < self.HEADER_SIZE:
            raise ClassSelectionError(
                "The result project is too small to contain a SIMPLE header."
            )
        count = self.MAX_SEGMENTS * self.HEADER_VALUES_PER_SEGMENT
        return struct.unpack("<" + ("q" * count), header)

    def segments(self):
        if self._segments is not None:
            return list(self._segments)

        segments = []
        project_size = os.path.getsize(self.project_path)
        values = self._header_values()
        for offset in range(0, len(values), self.HEADER_VALUES_PER_SEGMENT):
            index = (offset // self.HEADER_VALUES_PER_SEGMENT) + 1
            fromp, top, record_width, record_count, first_byte = values[
                offset : offset + self.HEADER_VALUES_PER_SEGMENT
            ]
            if record_width <= 0 or record_count <= 0:
                continue
            if fromp <= 0 or top < fromp or record_count != (top - fromp + 1):
                raise ClassSelectionError(
                    f"SIMPLE segment {index} has inconsistent record bounds."
                )
            if first_byte < self.HEADER_SIZE + 1:
                raise ClassSelectionError(
                    f"SIMPLE segment {index} has an invalid data offset."
                )
            if first_byte - 1 + (record_width * record_count) > project_size:
                raise ClassSelectionError(
                    f"SIMPLE segment {index} extends beyond the project file."
                )
            segments.append(
                SIMPLEProjectSegment(
                    index=index,
                    oritype=self.SEGMENT_NAMES.get(index, f"segment{index}"),
                    fromp=fromp,
                    top=top,
                    record_width=record_width,
                    record_count=record_count,
                    first_data_byte=first_byte,
                )
            )
        self._segments = segments
        return list(segments)

    def get_segment(self, oritype):
        return next(
            (segment for segment in self.segments() if segment.oritype == oritype),
            None,
        )

    def read_records(self, oritype):
        segment = self.get_segment(oritype)
        if segment is None:
            return []
        if segment.index in self.PARTICLE_SEGMENTS:
            raise ClassSelectionError(
                f"SIMPLE segment {oritype} is not a text metadata segment."
            )

        offset = segment.first_data_byte - 1
        records = []
        with open(self.project_path, "rb") as project_file:
            project_file.seek(offset)
            for _ in range(segment.record_count):
                raw_record = project_file.read(segment.record_width)
                if len(raw_record) != segment.record_width:
                    raise ClassSelectionError(
                        f"SIMPLE segment {oritype} ends unexpectedly."
                    )
                records.append(
                    self.parse_record(
                        raw_record.decode("utf-8", errors="replace").strip()
                    )
                )
        return records

    @classmethod
    def parse_record(cls, record):
        """Parse key/value metadata while preserving the first duplicate value."""
        parsed = {}
        for token in record.split():
            if "=" not in token:
                continue
            key, value = token.split("=", 1)
            if key not in parsed:
                parsed[key] = cls._coerce_value(value)
        return parsed

    @staticmethod
    def _coerce_value(value):
        if value == "":
            return value
        lower_value = value.lower()
        if lower_value in ("yes", "no"):
            return lower_value
        try:
            if any(character in value for character in (".", "e", "E")):
                return float(value)
            return int(value)
        except ValueError:
            return value


@dataclass(frozen=True)
class BatchClassSelection:
    """Validated server and browser data for one batch class-average stack."""

    stack_path: str
    stack_name: str
    classes: tuple
    initial_selected_class_ids: tuple
    storage_key: str
    stack_size: int
    stack_mtime_ns: int
    width: int
    height: int

    def browser_data(self):
        return {
            "classes": list(self.classes),
            "initial_selected_class_ids": list(self.initial_selected_class_ids),
            "selection_storage_key": self.storage_key,
        }


def _record_value(record, key):
    if not isinstance(record, dict):
        return None
    if key in record:
        return record[key]
    lower_key = key.lower()
    return next(
        (
            value
            for record_key, value in record.items()
            if str(record_key).lower() == lower_key
        ),
        None,
    )


def _numeric_record_value(record, key):
    value = _record_value(record, key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _display_number(value):
    if value is None:
        return None
    integer = int(value)
    return integer if value == integer else value


def _canonical_class_id(record, ordinal):
    value = _numeric_record_value(record, "class")
    if value is None:
        return ordinal + 1
    integer = int(value)
    if value != integer or integer <= 0:
        raise ClassSelectionError(
            f"cls2D record {ordinal + 1} has an invalid class identifier."
        )
    return integer


def _safe_regular_file(path, project_root):
    if not isinstance(path, (str, os.PathLike)) or not isinstance(
        project_root, (str, os.PathLike)
    ):
        return None
    try:
        resolved_path = Path(path).resolve(strict=True)
        resolved_root = Path(project_root).resolve(strict=True)
        contained = os.path.commonpath((resolved_path, resolved_root)) == str(
            resolved_root
        )
    except (OSError, ValueError):
        return None
    if not contained or not resolved_path.is_file():
        return None
    return resolved_path


def _class_average_stack(project_path, project_root, out_records):
    candidates = [
        record
        for record in out_records
        if str(_record_value(record, "imgkind") or "").lower() == "cavg"
    ]
    if len(candidates) != 1:
        raise ClassSelectionError(
            "The result project must contain exactly one class-average output stack."
        )

    record = candidates[0]
    raw_stack = _record_value(record, "stk")
    if not isinstance(raw_stack, str) or not raw_stack.strip():
        raise ClassSelectionError(
            "The class-average output does not identify an MRC stack."
        )

    project_directory = Path(project_path).parent
    raw_stack_path = Path(raw_stack.strip())
    paths = [
        raw_stack_path
        if raw_stack_path.is_absolute()
        else project_directory / raw_stack_path
    ]
    raw_stack_directory = _record_value(record, "stkpath")
    if isinstance(raw_stack_directory, str) and raw_stack_directory.strip():
        stack_directory = Path(raw_stack_directory.strip())
        if not stack_directory.is_absolute():
            stack_directory = project_directory / stack_directory
        stack_component = (
            raw_stack_path.name if raw_stack_path.is_absolute() else raw_stack_path
        )
        paths.append(stack_directory / stack_component)

    for path in paths:
        safe_path = _safe_regular_file(path, project_root)
        if safe_path is not None and safe_path.suffix.lower() in (".mrc", ".mrcs"):
            return safe_path
    raise ClassSelectionError(
        "The class-average MRC stack is missing or outside the selected project."
    )


@lru_cache(maxsize=64)
def _load_batch_class_selection(
    project_path,
    project_root,
    job_id,
    project_size,
    project_mtime_ns,
):
    """Build a dataset cached by the immutable result project's fingerprint."""
    del project_size, project_mtime_ns
    safe_project = _safe_regular_file(project_path, project_root)
    if safe_project is None or safe_project.suffix.lower() != ".simple":
        raise ClassSelectionError("The result project is unavailable or unsafe.")

    reader = SIMPLEProjectFileReader(safe_project)
    class_records = reader.read_records("cls2D")
    if not class_records:
        raise ClassSelectionError("The result project has no 2D class metadata.")
    stack_path = _class_average_stack(
        safe_project,
        project_root,
        reader.read_records("out"),
    )
    stack_info = read_mrc_stack_info(stack_path)
    if stack_info is None:
        raise ClassSelectionError("The class-average MRC stack is invalid.")
    if stack_info.count != len(class_records):
        raise ClassSelectionError(
            "The class-average stack and cls2D metadata contain different class counts."
        )

    classes = []
    selected_ids = []
    seen_ids = set()
    for ordinal, record in enumerate(class_records):
        class_id = _canonical_class_id(record, ordinal)
        if class_id in seen_ids:
            raise ClassSelectionError(
                f"The result project contains duplicate class ID {class_id}."
            )
        seen_ids.add(class_id)
        state = _numeric_record_value(record, "state")
        if state is None or state > 0:
            selected_ids.append(class_id)
        classes.append(
            {
                "class_id": class_id,
                "stack_index": ordinal + 1,
                "population": _display_number(
                    _numeric_record_value(record, "pop")
                ),
                "resolution": _numeric_record_value(record, "res"),
                "state": _display_number(state),
            }
        )

    project_stat = safe_project.stat()
    stack_stat = stack_path.stat()
    revision = hashlib.sha256(
        (
            f"{safe_project}:{project_stat.st_size}:{project_stat.st_mtime_ns}:"
            f"{stack_path}:{stack_stat.st_size}:{stack_stat.st_mtime_ns}"
        ).encode("utf-8")
    ).hexdigest()[:16]
    return BatchClassSelection(
        stack_path=str(stack_path),
        stack_name=stack_path.name,
        classes=tuple(classes),
        initial_selected_class_ids=tuple(sorted(selected_ids)),
        storage_key=f"nice.batch-class-selection.v1.{job_id}.{revision}",
        stack_size=stack_stat.st_size,
        stack_mtime_ns=stack_stat.st_mtime_ns,
        width=stack_info.width,
        height=stack_info.height,
    )


def load_batch_class_selection(project_path, project_root, job_id):
    """Return a validated class-selection dataset without changing source files."""
    safe_project = _safe_regular_file(project_path, project_root)
    if safe_project is None or safe_project.suffix.lower() != ".simple":
        raise ClassSelectionError("The result project is unavailable or unsafe.")
    project_stat = safe_project.stat()
    cache_args = (
        str(safe_project),
        str(Path(project_root).resolve(strict=True)),
        job_id,
        project_stat.st_size,
        project_stat.st_mtime_ns,
    )
    selection = _load_batch_class_selection(*cache_args)
    try:
        stack_stat = Path(selection.stack_path).stat()
    except OSError:
        _load_batch_class_selection.cache_clear()
        raise ClassSelectionError("The class-average MRC stack is unavailable.")
    if (
        stack_stat.st_size != selection.stack_size
        or stack_stat.st_mtime_ns != selection.stack_mtime_ns
    ):
        _load_batch_class_selection.cache_clear()
        selection = _load_batch_class_selection(*cache_args)
    return selection


def batch_class_selection_available(project_path, project_root):
    """Check for non-empty cls2D metadata without loading pixels or writing files."""
    safe_project = _safe_regular_file(project_path, project_root)
    if safe_project is None or safe_project.suffix.lower() != ".simple":
        return False
    try:
        segment = SIMPLEProjectFileReader(safe_project).get_segment("cls2D")
    except (ClassSelectionError, OSError, OverflowError, struct.error):
        return False
    return segment is not None and segment.record_count > 0


def deselected_class_ids(selection, selected_class_ids):
    """Validate browser state and return canonical deselected class IDs."""
    if not isinstance(selected_class_ids, list):
        raise ClassSelectionError("Selection data is missing or invalid.")
    if any(
        isinstance(value, bool) or not isinstance(value, int)
        for value in selected_class_ids
    ):
        raise ClassSelectionError("Selection contains a non-integer class ID.")

    known_ids = {entry["class_id"] for entry in selection.classes}
    selected_ids = set(selected_class_ids)
    unknown_ids = selected_ids - known_ids
    if unknown_ids:
        raise ClassSelectionError(
            f"Selection contains unknown class IDs: {sorted(unknown_ids)}."
        )
    return sorted(known_ids - selected_ids)
