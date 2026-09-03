import os
import struct
import tempfile
from unittest import mock

from django.test import SimpleTestCase

from ..data_structures import mrc as mrc_helpers
from ..data_structures.mrc import (
    read_mrc_stack_info,
    render_mrc_movie_png,
    render_mrc_movie_webp,
    render_mrc_particle_png,
)


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

    def test_reads_stack_layout_without_creating_thumbnail_files(self):
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

        with mock.patch.object(
            mrc_helpers.mrcfile,
            "mmap",
            wraps=mrc_helpers.mrcfile.mmap,
        ) as mrc_mmap, mock.patch.object(
            mrc_helpers.Image,
            "fromarray",
            wraps=mrc_helpers.Image.fromarray,
        ) as image_fromarray:
            thumbnail = render_mrc_particle_png(stack_path, 2)
        mrc_mmap.assert_called_once_with(stack_path, mode="r", permissive=True)
        image_fromarray.assert_called_once()

        self.assertEqual((info.width, info.height, info.count, info.mode), (2, 2, 2, 2))
        self.assertEqual(info.data_offset, 1040)
        self.assertEqual(thumbnail[:8], b"\x89PNG\r\n\x1a\n")
        self.assertEqual(struct.unpack(">2I", thumbnail[16:24]), (2, 2))
        self.assertEqual(set(os.listdir(self.tempdir.name)), before)

    def test_renders_float16_extract_stack_and_rejects_invalid_index(self):
        stack_path = os.path.join(self.tempdir.name, "particles_float16.mrcs")
        _write_mrc_stack(
            stack_path,
            width=2,
            height=2,
            images=((0.0, 0.25, 0.5, 1.0),),
            mode=12,
        )

        info = read_mrc_stack_info(stack_path)

        self.assertEqual(info.mode, 12)
        self.assertIsNotNone(render_mrc_particle_png(stack_path, 1))
        self.assertIsNone(render_mrc_particle_png(stack_path, 0))
        self.assertIsNone(render_mrc_particle_png(stack_path, 2))

    def test_renders_sampled_movie_average_without_creating_a_thumbnail(self):
        movie_path = os.path.join(self.tempdir.name, "movie.mrcs")
        _write_mrc_stack(
            movie_path,
            width=4,
            height=2,
            images=(
                (0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0),
                (0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
            ),
        )

        before = set(os.listdir(self.tempdir.name))
        with mock.patch.object(
            mrc_helpers.mrcfile,
            "mmap",
            wraps=mrc_helpers.mrcfile.mmap,
        ) as mrc_mmap:
            thumbnail = render_mrc_movie_png(movie_path, max_size=2)

        mrc_mmap.assert_called_once_with(movie_path, mode="r", permissive=True)
        self.assertEqual(thumbnail[:8], b"\x89PNG\r\n\x1a\n")
        self.assertEqual(struct.unpack(">2I", thumbnail[16:24]), (2, 1))
        webp_thumbnail = render_mrc_movie_webp(movie_path, max_size=2)
        self.assertEqual(webp_thumbnail[:4], b"RIFF")
        self.assertEqual(webp_thumbnail[8:12], b"WEBP")
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
        self.assertIsNone(render_mrc_particle_png(stack_path, 1))
