import os
import struct
import tempfile
from unittest import mock

import numpy as np
from django.test import SimpleTestCase

from ..data_structures import mrc as mrc_helpers
from ..data_structures.mrc import read_mrc_stack_info, read_mrc_volume_info


def _write_mrc_stack(path, width, height, images, mode=2, extended_header_size=0):
    value_format = {2: "f", 12: "e"}[mode]
    header = bytearray(1024)
    struct.pack_into("<4i", header, 0, width, height, len(images), mode)
    struct.pack_into("<i", header, 92, extended_header_size)
    header[208:212] = b"MAP "
    header[212:216] = b"DA\x00\x00"
    with open(path, "wb") as stack_file:
        stack_file.write(header)
        stack_file.write(b"\x00" * extended_header_size)
        for image in images:
            stack_file.write(struct.pack(f"<{len(image)}{value_format}", *image))


class MRCStackTests(SimpleTestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)

    def test_reads_stack_layout_without_creating_files(self):
        stack_path = os.path.join(self.tempdir.name, "particles.mrcs")
        _write_mrc_stack(
            stack_path,
            width=2,
            height=2,
            images=((0.0, 1.0, 2.0, 3.0), (4.0, 5.0, 6.0, 7.0)),
            extended_header_size=16,
        )

        before = set(os.listdir(self.tempdir.name))
        with mock.patch.object(
            mrc_helpers.mrcfile,
            "open",
            wraps=mrc_helpers.mrcfile.open,
        ) as mrc_open:
            info = read_mrc_stack_info(stack_path)

        mrc_open.assert_called_once_with(
            stack_path,
            mode="r",
            permissive=True,
            header_only=True,
        )
        self.assertEqual((info.width, info.height, info.count, info.mode), (2, 2, 2, 2))
        self.assertEqual(info.data_offset, 1040)
        self.assertEqual(set(os.listdir(self.tempdir.name)), before)

    def test_rejects_truncated_stack(self):
        stack_path = os.path.join(self.tempdir.name, "truncated.mrcs")
        _write_mrc_stack(
            stack_path,
            width=2,
            height=2,
            images=((0.0, 1.0, 2.0, 3.0),),
        )
        with open(stack_path, "rb+") as stack_file:
            stack_file.truncate(1024)

        self.assertIsNone(read_mrc_stack_info(stack_path))

    def test_reads_volume_metadata(self):
        volume_path = os.path.join(self.tempdir.name, "volume.mrc")
        source = np.arange(10 * 9 * 8, dtype=np.float32).reshape((10, 9, 8))
        with mrc_helpers.mrcfile.new(volume_path) as volume:
            volume.set_data(source)
            volume.voxel_size = (1.5, 2.0, 2.5)
            volume.update_header_stats()

        info = read_mrc_volume_info(volume_path)

        self.assertEqual((info.width, info.height, info.depth), (8, 9, 10))
        self.assertEqual(info.voxel_size, (1.5, 2.0, 2.5))

    def test_volume_metadata_rejects_a_single_image_stack(self):
        stack_path = os.path.join(self.tempdir.name, "single.mrc")
        _write_mrc_stack(
            stack_path,
            width=2,
            height=2,
            images=((0.0, 1.0, 2.0, 3.0),),
        )

        self.assertIsNone(read_mrc_volume_info(stack_path))
