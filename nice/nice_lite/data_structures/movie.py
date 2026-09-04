"""On-demand movie thumbnail rendering for NICE batch imports."""

import os
import struct

import numpy as np
from PIL import Image

from .mrc import (
    MAX_MOVIE_PREVIEW_FRAMES,
    MAX_MRC_DIMENSION,
    read_mrc_stack_info,
    render_grayscale_array_png,
    render_grayscale_array_webp,
    render_mrc_movie_png,
    render_mrc_movie_webp,
)


MRC_MOVIE_EXTENSIONS = frozenset((".mrc", ".mrcs"))
TIFF_MOVIE_EXTENSIONS = frozenset((".tif", ".tiff"))
TIFF_DIMENSION_EXTENSIONS = TIFF_MOVIE_EXTENSIONS | frozenset((".eer",))
SUPPORTED_MOVIE_EXTENSIONS = MRC_MOVIE_EXTENSIONS | TIFF_MOVIE_EXTENSIONS
_MAX_TIFF_IFD_ENTRIES = 4096
_TIFF_INTEGER_FORMATS = {
    3: "H",  # SHORT
    4: "I",  # LONG
    16: "Q",  # LONG8
}


def _read_tiff_dimensions(path):
    """Read width and height from the first classic TIFF or BigTIFF IFD."""
    try:
        with open(path, "rb") as movie:
            header = movie.read(16)
            if header[:2] == b"II":
                byte_order = "<"
            elif header[:2] == b"MM":
                byte_order = ">"
            else:
                return None

            version = struct.unpack(byte_order + "H", header[2:4])[0]
            if version == 42:
                if len(header) < 8:
                    return None
                ifd_offset = struct.unpack(byte_order + "I", header[4:8])[0]
                count_format = "H"
                count_size = 2
                entry_size = 12
                offset_size = 4
            elif version == 43:
                if len(header) < 16:
                    return None
                stored_offset_size, reserved = struct.unpack(
                    byte_order + "HH",
                    header[4:8],
                )
                if stored_offset_size != 8 or reserved != 0:
                    return None
                ifd_offset = struct.unpack(byte_order + "Q", header[8:16])[0]
                count_format = "Q"
                count_size = 8
                entry_size = 20
                offset_size = 8
            else:
                return None

            file_size = os.fstat(movie.fileno()).st_size
            if ifd_offset < 8 or ifd_offset + count_size > file_size:
                return None
            movie.seek(ifd_offset)
            entry_count_data = movie.read(count_size)
            if len(entry_count_data) != count_size:
                return None
            entry_count = struct.unpack(
                byte_order + count_format,
                entry_count_data,
            )[0]
            if not 0 < entry_count <= _MAX_TIFF_IFD_ENTRIES:
                return None
            if (
                ifd_offset
                + count_size
                + (entry_count * entry_size)
                + offset_size
                > file_size
            ):
                return None

            dimensions = {}
            for _ in range(entry_count):
                entry = movie.read(entry_size)
                if len(entry) != entry_size:
                    return None
                if version == 42:
                    tag, value_type, value_count = struct.unpack(
                        byte_order + "HHI",
                        entry[:8],
                    )
                    value_bytes = entry[8:12]
                else:
                    tag, value_type, value_count = struct.unpack(
                        byte_order + "HHQ",
                        entry[:12],
                    )
                    value_bytes = entry[12:20]

                value_format = _TIFF_INTEGER_FORMATS.get(value_type)
                if tag not in (256, 257) or value_count != 1 or value_format is None:
                    continue
                value_size = struct.calcsize(value_format)
                value = struct.unpack(
                    byte_order + value_format,
                    value_bytes[:value_size],
                )[0]
                dimensions[tag] = int(value)

            width = dimensions.get(256)
            height = dimensions.get(257)
    except (OSError, OverflowError, struct.error, TypeError, ValueError):
        return None

    if (
        not isinstance(width, int)
        or not isinstance(height, int)
        or not 0 < width <= MAX_MRC_DIMENSION
        or not 0 < height <= MAX_MRC_DIMENSION
    ):
        return None
    return width, height


def read_movie_dimensions(path):
    """Return source-frame width and height without reading movie pixels."""
    if not isinstance(path, str):
        return None

    extension = os.path.splitext(path)[1].lower()
    if extension in MRC_MOVIE_EXTENSIONS:
        info = read_mrc_stack_info(path)
        return (info.width, info.height) if info is not None else None
    if extension in TIFF_DIMENSION_EXTENSIONS:
        return _read_tiff_dimensions(path)
    return None


def movie_preview_supported(path):
    """Return whether the movie extension has an on-demand renderer."""
    return (
        isinstance(path, str)
        and os.path.splitext(path)[1].lower() in SUPPORTED_MOVIE_EXTENSIONS
    )


def _read_tiff_movie_preview(path, max_size):
    """Return an averaged, sampled TIFF movie preview array."""
    try:
        with Image.open(path) as movie:
            width, height = movie.size
            frame_count = int(getattr(movie, "n_frames", 1))
            if (
                not 0 < width <= MAX_MRC_DIMENSION
                or not 0 < height <= MAX_MRC_DIMENSION
                or frame_count < 1
            ):
                return None

            scale = min(1.0, max_size / max(width, height))
            output_size = (
                max(1, round(width * scale)),
                max(1, round(height * scale)),
            )
            frame_indices = np.unique(np.linspace(
                0,
                frame_count - 1,
                num=min(frame_count, MAX_MOVIE_PREVIEW_FRAMES),
                dtype=np.intp,
            ))
            averaged_values = np.zeros(
                (output_size[1], output_size[0]),
                dtype=np.float64,
            )
            for frame_index in frame_indices:
                movie.seek(int(frame_index))
                # Convert 16-bit detector counts before resizing. Resizing the
                # native I;16 image first quantizes this data to an almost
                # binary image and produces a visually blank thumbnail.
                frame = movie.convert("F")
                if frame.size != output_size:
                    frame = frame.resize(output_size, Image.Resampling.BOX)
                averaged_values += np.asarray(frame, dtype=np.float64)
            averaged_values /= frame_indices.size
    except (EOFError, OSError, OverflowError, TypeError, ValueError):
        return None

    return averaged_values


def render_movie_png(path, max_size=160):
    """Render a supported movie on demand without creating a thumbnail file."""
    if (
        not isinstance(path, str)
        or not isinstance(max_size, int)
        or isinstance(max_size, bool)
        or max_size < 1
    ):
        return None

    extension = os.path.splitext(path)[1].lower()
    if extension in MRC_MOVIE_EXTENSIONS:
        return render_mrc_movie_png(path, max_size=max_size)
    if extension in TIFF_MOVIE_EXTENSIONS:
        values = _read_tiff_movie_preview(path, max_size)
        return render_grayscale_array_png(values) if values is not None else None
    return None


def render_movie_webp(path, max_size=160):
    """Render a supported movie as an in-memory WebP thumbnail."""
    if (
        not isinstance(path, str)
        or not isinstance(max_size, int)
        or isinstance(max_size, bool)
        or max_size < 1
    ):
        return None

    extension = os.path.splitext(path)[1].lower()
    if extension in MRC_MOVIE_EXTENSIONS:
        return render_mrc_movie_webp(path, max_size=max_size)
    if extension in TIFF_MOVIE_EXTENSIONS:
        values = _read_tiff_movie_preview(path, max_size)
        return render_grayscale_array_webp(values) if values is not None else None
    return None
