import os
import struct
import tempfile
from unittest import mock

import numpy as np
from django.test import SimpleTestCase
from PIL import Image

from ..data_structures.movie import (
    movie_preview_supported,
    read_movie_dimensions,
    render_movie_webp,
)


class MovieThumbnailTests(SimpleTestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)

    def test_renders_multipage_tiff_without_creating_a_thumbnail_file(self):
        movie_path = os.path.join(self.tempdir.name, "movie.tiff")
        frames = [
            Image.fromarray(np.array([[0, 1], [2, 3]], dtype=np.uint16)),
            Image.fromarray(np.array([[3, 2], [1, 0]], dtype=np.uint16)),
        ]
        frames[0].save(movie_path, save_all=True, append_images=frames[1:])

        before = set(os.listdir(self.tempdir.name))
        thumbnail = render_movie_webp(movie_path)

        self.assertEqual(thumbnail[:4], b"RIFF")
        self.assertEqual(thumbnail[8:12], b"WEBP")
        self.assertEqual(set(os.listdir(self.tempdir.name)), before)

    def test_converts_sixteen_bit_tiff_data_before_downsampling(self):
        movie_path = os.path.join(self.tempdir.name, "detector_counts.tiff")
        Image.fromarray(np.array([
            [0, 1, 0, 1],
            [1, 1, 1, 1],
            [0, 1, 25, 1],
            [1, 1, 1, 1],
        ], dtype=np.uint16)).save(movie_path)
        operations = []
        original_convert = Image.Image.convert
        original_resize = Image.Image.resize

        def record_convert(image, *args, **kwargs):
            operations.append("convert")
            return original_convert(image, *args, **kwargs)

        def record_resize(image, *args, **kwargs):
            operations.append("resize")
            return original_resize(image, *args, **kwargs)

        with (
            mock.patch.object(Image.Image, "convert", record_convert),
            mock.patch.object(Image.Image, "resize", record_resize),
        ):
            thumbnail = render_movie_webp(movie_path, max_size=2)

        self.assertIsNotNone(thumbnail)
        self.assertLess(operations.index("convert"), operations.index("resize"))

    def test_reads_tiff_dimensions_without_loading_movie_pixels(self):
        movie_path = os.path.join(self.tempdir.name, "movie.tiff")
        Image.fromarray(np.zeros((3, 5), dtype=np.uint16)).save(movie_path)

        with mock.patch.object(
            Image.Image,
            "load",
            side_effect=AssertionError("movie pixels should remain lazy"),
        ):
            dimensions = read_movie_dimensions(movie_path)

        self.assertEqual(dimensions, (5, 3))

    def test_reads_big_tiff_eer_dimensions(self):
        movie_path = os.path.join(self.tempdir.name, "movie.eer")
        with open(movie_path, "wb") as movie:
            movie.write(struct.pack("<2sHHHQ", b"II", 43, 8, 0, 16))
            movie.write(struct.pack("<Q", 2))
            movie.write(struct.pack("<HHQQ", 256, 4, 1, 4096))
            movie.write(struct.pack("<HHQQ", 257, 4, 1, 4096))
            movie.write(struct.pack("<Q", 0))

        self.assertEqual(read_movie_dimensions(movie_path), (4096, 4096))

    def test_rejects_unsupported_movie_formats(self):
        self.assertFalse(movie_preview_supported("movie.eer"))
        self.assertIsNone(render_movie_webp("movie.eer"))
