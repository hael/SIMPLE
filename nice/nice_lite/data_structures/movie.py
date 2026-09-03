"""On-demand movie thumbnail rendering for NICE batch imports."""

import os

import numpy as np
from PIL import Image

from .mrc import (
    MAX_MOVIE_PREVIEW_FRAMES,
    MAX_MRC_DIMENSION,
    render_grayscale_array_png,
    render_grayscale_array_webp,
    render_mrc_movie_png,
    render_mrc_movie_webp,
)


MRC_MOVIE_EXTENSIONS = frozenset((".mrc", ".mrcs"))
TIFF_MOVIE_EXTENSIONS = frozenset((".tif", ".tiff"))
SUPPORTED_MOVIE_EXTENSIONS = MRC_MOVIE_EXTENSIONS | TIFF_MOVIE_EXTENSIONS


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
