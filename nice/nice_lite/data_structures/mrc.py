"""Validated MRC stack and volume metadata for NICE."""

import os
import warnings
from dataclasses import dataclass
import mrcfile
import numpy as np


MRC_HEADER_SIZE = 1024
MAX_MRC_DIMENSION = 32768
MAX_MRC_IMAGES = 100000000
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
