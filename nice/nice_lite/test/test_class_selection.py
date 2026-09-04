import os
import struct
import tempfile

from django.test import SimpleTestCase

from ..data_structures.class_selection import (
    ClassSelectionError,
    SIMPLEProjectFileReader,
    deselected_class_ids,
    load_batch_class_selection,
)
from .test_mrc import _write_mrc_stack


def _write_project(path, class_records, out_records):
    segment_records = {4: class_records, 7: out_records}
    header = [0] * (
        SIMPLEProjectFileReader.MAX_SEGMENTS
        * SIMPLEProjectFileReader.HEADER_VALUES_PER_SEGMENT
    )
    payloads = []
    first_data_byte = SIMPLEProjectFileReader.HEADER_SIZE + 1
    for segment_index in sorted(segment_records):
        records = [record.encode("utf-8") for record in segment_records[segment_index]]
        record_width = max(len(record) for record in records)
        offset = (segment_index - 1) * SIMPLEProjectFileReader.HEADER_VALUES_PER_SEGMENT
        header[offset : offset + 5] = [
            1,
            len(records),
            record_width,
            len(records),
            first_data_byte,
        ]
        payload = b"".join(record.ljust(record_width, b" ") for record in records)
        payloads.append(payload)
        first_data_byte += len(payload)

    with open(path, "wb") as project_file:
        project_file.write(struct.pack("<" + ("q" * len(header)), *header))
        for payload in payloads:
            project_file.write(payload)


class BatchClassSelectionTests(SimpleTestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.stack_path = os.path.join(self.tempdir.name, "classes.mrcs")
        self.project_path = os.path.join(self.tempdir.name, "workspace.simple")
        _write_mrc_stack(
            self.stack_path,
            width=2,
            height=2,
            images=(
                (0.0, 1.0, 2.0, 3.0),
                (4.0, 5.0, 6.0, 7.0),
                (8.0, 9.0, 10.0, 11.0),
            ),
        )

    def _write_selection_project(self, stack_path=None):
        _write_project(
            self.project_path,
            class_records=[
                "class=1 pop=10 res=4.5 state=1",
                "class=2 pop=20 res=7.0 state=0",
                "class=3 pop=30 res=5.5 state=1",
            ],
            out_records=[f"imgkind=cavg stk={stack_path or self.stack_path}"],
        )

    def test_loads_all_classes_and_initializes_selection_from_project_state(self):
        self._write_selection_project()

        selection = load_batch_class_selection(
            self.project_path,
            self.tempdir.name,
            job_id=7,
        )

        self.assertEqual(selection.stack_path, self.stack_path)
        self.assertEqual(selection.initial_selected_class_ids, (1, 3))
        self.assertEqual(
            [entry["class_id"] for entry in selection.classes],
            [1, 2, 3],
        )
        self.assertEqual(selection.classes[1]["population"], 20)
        self.assertEqual((selection.width, selection.height), (2, 2))
        self.assertIn("nice.batch-class-selection.v1.7", selection.storage_key)

    def test_rejects_class_stack_outside_selected_project(self):
        with tempfile.TemporaryDirectory() as outside_directory:
            outside_stack = os.path.join(outside_directory, "outside.mrcs")
            _write_mrc_stack(
                outside_stack,
                width=2,
                height=2,
                images=((0.0, 1.0, 2.0, 3.0),) * 3,
            )
            self._write_selection_project(stack_path=outside_stack)

            with self.assertRaisesRegex(ClassSelectionError, "outside"):
                load_batch_class_selection(
                    self.project_path,
                    self.tempdir.name,
                    job_id=7,
                )

    def test_deselection_export_uses_canonical_one_based_ids(self):
        self._write_selection_project()
        selection = load_batch_class_selection(
            self.project_path,
            self.tempdir.name,
            job_id=7,
        )

        self.assertEqual(deselected_class_ids(selection, [3, 1]), [2])
        with self.assertRaisesRegex(ClassSelectionError, "unknown class IDs"):
            deselected_class_ids(selection, [1, 4])

    def test_parser_preserves_first_duplicate_value(self):
        record = SIMPLEProjectFileReader.parse_record("stk=classes.mrcs stk=0")

        self.assertEqual(record["stk"], "classes.mrcs")
