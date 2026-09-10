"""Lazy MRC stack metadata and particle PNG rendering for NICE."""

import os
import warnings
from dataclasses import dataclass
from io import BytesIO

import mrcfile
import numpy as np
from PIL import Image


MRC_HEADER_SIZE = 1024
MAX_MRC_DIMENSION = 32768
MAX_MRC_IMAGES = 100000000
MAX_RENDER_PIXELS = 4194304
MAX_MOVIE_PREVIEW_FRAMES = 8
SUPPORTED_MRC_MODES = frozenset((0, 1, 2, 6, 12))


@dataclass(frozen=True)
class MRCStackInfo:
    """Validated layout information needed to address one MRC stack image."""

    width: int
    height: int
    count: int
    mode: int
    data_offset: int


@dataclass(frozen=True)
class MRCVolumeInfo:
    """Validated metadata for one three-dimensional MRC density map."""

    width: int
    height: int
    depth: int
    mode: int
    data_offset: int
    voxel_size: tuple
    minimum: float
    maximum: float


@dataclass(frozen=True)
class MRCVolumePayload:
    """A bounded, normalized volume ready for a browser 3D texture."""

    data: bytes
    width: int
    height: int
    depth: int
    source_width: int
    source_height: int
    source_depth: int
    voxel_size: tuple
    minimum: float
    maximum: float


def read_mrc_stack_info(path):
    """Return validated MRC stack metadata without reading particle pixels."""
    try:
        file_size = os.path.getsize(path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with mrcfile.open(
                path,
                mode="r",
                permissive=True,
                header_only=True,
            ) as stack:
                width = int(stack.header.nx)
                height = int(stack.header.ny)
                count = int(stack.header.nz)
                mode = int(stack.header.mode)
                extended_header_size = int(stack.header.nsymbt)
                dtype = mrcfile.utils.data_dtype_from_header(stack.header)
    except (OSError, OverflowError, TypeError, ValueError):
        return None

    if (
        mode not in SUPPORTED_MRC_MODES
        or not 0 < width <= MAX_MRC_DIMENSION
        or not 0 < height <= MAX_MRC_DIMENSION
        or not 0 < count <= MAX_MRC_IMAGES
        or extended_header_size < 0
    ):
        return None

    data_offset = MRC_HEADER_SIZE + extended_header_size
    expected_size = data_offset + (width * height * count * dtype.itemsize)
    if data_offset > file_size or expected_size > file_size:
        return None

    return MRCStackInfo(
        width=width,
        height=height,
        count=count,
        mode=mode,
        data_offset=data_offset,
    )


def read_mrc_volume_info(path):
    """Return validated 3D-map metadata without reading density voxels."""
    stack_info = read_mrc_stack_info(path)
    if stack_info is None or stack_info.count <= 1:
        return None

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with mrcfile.open(
                path,
                mode="r",
                permissive=True,
                header_only=True,
            ) as volume:
                voxel_size = tuple(
                    float(axis)
                    for axis in (
                        volume.voxel_size.x,
                        volume.voxel_size.y,
                        volume.voxel_size.z,
                    )
                )
                minimum = float(volume.header.dmin)
                maximum = float(volume.header.dmax)
    except (OSError, OverflowError, TypeError, ValueError):
        return None

    if not all(np.isfinite(axis) and axis > 0.0 for axis in voxel_size):
        voxel_size = (1.0, 1.0, 1.0)
    if not np.isfinite(minimum) or not np.isfinite(maximum):
        minimum = maximum = 0.0

    return MRCVolumeInfo(
        width=stack_info.width,
        height=stack_info.height,
        depth=stack_info.count,
        mode=stack_info.mode,
        data_offset=stack_info.data_offset,
        voxel_size=voxel_size,
        minimum=minimum,
        maximum=maximum,
    )


def _sample_volume_axis(length, max_dimension):
    output_length = min(length, max_dimension)
    return np.minimum(
        length - 1,
        ((np.arange(output_length) + 0.5) * length / output_length).astype(np.intp),
    )


def read_mrc_volume_payload(path, max_dimension=128):
    """Read a bounded 3D texture, normalized to the source density range."""
    if (
        not isinstance(max_dimension, int)
        or isinstance(max_dimension, bool)
        or not 8 <= max_dimension <= 256
    ):
        return None

    info = read_mrc_volume_info(path)
    if info is None:
        return None

    source_slices = _sample_volume_axis(info.depth, max_dimension)
    source_rows = _sample_volume_axis(info.height, max_dimension)
    source_columns = _sample_volume_axis(info.width, max_dimension)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with mrcfile.mmap(path, mode="r", permissive=True) as volume:
                values = volume.data
                if (
                    values is None
                    or values.ndim != 3
                    or values.shape != (info.depth, info.height, info.width)
                ):
                    return None
                sampled = np.asarray(
                    values[np.ix_(source_slices, source_rows, source_columns)],
                    dtype=np.float32,
                )
    except (IndexError, OSError, OverflowError, TypeError, ValueError):
        return None

    finite_mask = np.isfinite(sampled)
    if not finite_mask.any():
        return None

    minimum = info.minimum
    maximum = info.maximum
    if maximum <= minimum:
        finite_values = sampled[finite_mask]
        minimum = float(np.min(finite_values))
        maximum = float(np.max(finite_values))

    dynamic_range = maximum - minimum
    if dynamic_range <= 0.0:
        texture = np.zeros(sampled.shape, dtype=np.uint8)
    else:
        safe_values = np.where(finite_mask, sampled, minimum)
        texture = np.rint(
            np.clip(255.0 * (safe_values - minimum) / dynamic_range, 0.0, 255.0)
        ).astype(np.uint8)

    return MRCVolumePayload(
        data=np.ascontiguousarray(texture).tobytes(order="C"),
        width=texture.shape[2],
        height=texture.shape[1],
        depth=texture.shape[0],
        source_width=info.width,
        source_height=info.height,
        source_depth=info.depth,
        voxel_size=info.voxel_size,
        minimum=minimum,
        maximum=maximum,
    )


def _encode_grayscale_image(pixels, image_format):
    """Encode a two-dimensional uint8 array as an in-memory image."""
    if (
        pixels.ndim != 2
        or pixels.dtype != np.uint8
        or image_format not in ("PNG", "WEBP")
    ):
        return None
    output = BytesIO()
    save_options = {"compress_level": 6} if image_format == "PNG" else {
        "quality": 80,
        "method": 4,
    }
    Image.fromarray(pixels).save(output, format=image_format, **save_options)
    return output.getvalue()


def _render_grayscale_array(values, image_format):
    """Normalize a two-dimensional numeric array and encode it."""
    if not isinstance(values, np.ndarray) or values.ndim != 2 or values.size == 0:
        return None

    finite_mask = np.isfinite(values)
    if finite_mask.any():
        finite_values = values[finite_mask].astype(np.float64, copy=False)
        minimum = float(np.min(finite_values))
        maximum = float(np.max(finite_values))
        mean = float(np.mean(finite_values))
        standard_deviation = float(np.std(finite_values))
        lower = max(minimum, mean - (3.0 * standard_deviation))
        upper = min(maximum, mean + (3.0 * standard_deviation))
        if upper <= lower:
            lower, upper = minimum, maximum
    else:
        lower = upper = 0.0

    dynamic_range = upper - lower
    if dynamic_range <= 0.0:
        pixels = np.full(values.shape, 127, dtype=np.uint8)
    else:
        safe_values = np.where(finite_mask, values, lower)
        pixels = np.rint(
            np.clip(255.0 * (safe_values - lower) / dynamic_range, 0.0, 255.0)
        ).astype(np.uint8)
    pixels[~finite_mask] = 0
    return _encode_grayscale_image(pixels, image_format)


def render_grayscale_array_png(values):
    """Normalize a two-dimensional numeric array and encode it as PNG."""
    return _render_grayscale_array(values, "PNG")


def render_grayscale_array_webp(values):
    """Normalize a two-dimensional numeric array and encode it as WebP."""
    return _render_grayscale_array(values, "WEBP")


def _sample_indices(width, height, max_size):
    """Return output dimensions and nearest-neighbour source indices."""
    scale = min(1.0, max_size / max(width, height))
    output_width = max(1, round(width * scale))
    output_height = max(1, round(height * scale))
    source_columns = np.minimum(
        width - 1,
        ((np.arange(output_width) + 0.5) * width / output_width).astype(np.intp),
    )
    source_rows = np.minimum(
        height - 1,
        ((np.arange(output_height) + 0.5) * height / output_height).astype(np.intp),
    )
    return source_rows, source_columns


def _read_particle(path, info, particle_index):
    """Copy one particle from an on-demand memory-mapped MRC data slice."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with mrcfile.mmap(path, mode="r", permissive=True) as stack:
                data = stack.data
                if data is None:
                    return None
                if data.ndim == 2:
                    if info.count != 1 or particle_index != 1:
                        return None
                    particle = data
                elif data.ndim == 3:
                    if data.shape != (info.count, info.height, info.width):
                        return None
                    particle = data[particle_index - 1]
                else:
                    return None
                if particle.shape != (info.height, info.width):
                    return None
                return np.array(particle, dtype=np.float32, copy=True)
    except (IndexError, OSError, OverflowError, TypeError, ValueError):
        return None


def render_mrc_particle_png(path, particle_index, max_size=160):
    """Render one 1-based stack image to PNG bytes without creating a file."""
    if (
        not isinstance(particle_index, int)
        or isinstance(particle_index, bool)
        or particle_index < 1
        or not isinstance(max_size, int)
        or isinstance(max_size, bool)
        or max_size < 1
    ):
        return None

    info = read_mrc_stack_info(path)
    if info is None or particle_index > info.count:
        return None

    pixel_count = info.width * info.height
    if pixel_count > MAX_RENDER_PIXELS:
        return None
    values = _read_particle(path, info, particle_index)
    if values is None or values.size != pixel_count:
        return None

    source_rows, source_columns = _sample_indices(
        info.width,
        info.height,
        max_size,
    )
    sampled_values = values[np.ix_(source_rows, source_columns)]
    return render_grayscale_array_png(sampled_values)


def _read_mrc_movie_preview(path, max_size):
    """Return an averaged, sampled MRC movie preview array."""
    if (
        not isinstance(max_size, int)
        or isinstance(max_size, bool)
        or max_size < 1
    ):
        return None

    info = read_mrc_stack_info(path)
    if info is None:
        return None

    source_rows, source_columns = _sample_indices(
        info.width,
        info.height,
        max_size,
    )
    frame_indices = np.unique(np.linspace(
        0,
        info.count - 1,
        num=min(info.count, MAX_MOVIE_PREVIEW_FRAMES),
        dtype=np.intp,
    ))
    averaged_values = np.zeros(
        (source_rows.size, source_columns.size),
        dtype=np.float64,
    )

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with mrcfile.mmap(path, mode="r", permissive=True) as movie:
                data = movie.data
                if data is None:
                    return None
                if data.ndim == 2:
                    if info.count != 1 or data.shape != (info.height, info.width):
                        return None
                    averaged_values += data[np.ix_(source_rows, source_columns)]
                elif data.ndim == 3:
                    if data.shape != (info.count, info.height, info.width):
                        return None
                    for frame_index in frame_indices:
                        averaged_values += data[frame_index][
                            np.ix_(source_rows, source_columns)
                        ]
                    averaged_values /= frame_indices.size
                else:
                    return None
    except (IndexError, OSError, OverflowError, TypeError, ValueError):
        return None

    return averaged_values


def render_mrc_movie_png(path, max_size=160):
    """Render an averaged MRC movie preview as an in-memory PNG."""
    values = _read_mrc_movie_preview(path, max_size)
    return render_grayscale_array_png(values) if values is not None else None


def render_mrc_movie_webp(path, max_size=160):
    """Render an averaged MRC movie preview as an in-memory WebP image."""
    values = _read_mrc_movie_preview(path, max_size)
    return render_grayscale_array_webp(values) if values is not None else None
