"""Prepare compact MRC volumes for the NICE Mol* viewer."""

import tempfile
from pathlib import Path

import mrcfile
import numpy as np
from django.http import HttpResponse
from django.utils.cache import patch_cache_control


VOLUME_BROWSER_CACHE_SECONDS = 86400


def volume_cache_version(path):
    """Identify the current source file for a versioned browser URL."""
    try:
        stat = Path(path).stat()
    except OSError:
        return None
    return f'{stat.st_mtime_ns}-{stat.st_size}'


def to_int8(data, minimum, maximum):
    values = np.abs(data) if np.iscomplexobj(data) else data
    values = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(values)
    result = np.zeros(values.shape, dtype=np.int8)
    if not finite.any():
        return result

    if not np.isfinite(minimum) or not np.isfinite(maximum) or maximum <= minimum:
        minimum, maximum = values[finite].min(), values[finite].max()
    if maximum > minimum:
        scaled = np.clip((values[finite] - minimum) / (maximum - minimum), 0, 1)
        result[finite] = np.rint(scaled * 255 - 128).astype(np.int8)
    return result


def spherical_mask(shape):
    center = (np.array(shape) - 1) / 2
    radius = min(shape) / 2
    z, y, x = np.ogrid[:shape[0], :shape[1], :shape[2]]
    return ((z - center[0]) ** 2
            + (y - center[1]) ** 2
            + (x - center[2]) ** 2 <= radius ** 2)


def compressed_volume_response(input_file, requested_version=None):
    """Convert a map to a gzip MRC response without keeping a copy."""
    source_version = volume_cache_version(input_file)
    with mrcfile.open(input_file) as source:
        mask = spherical_mask(source.data.shape)
        output_data = np.full(source.data.shape, -128, dtype=np.int8)
        # Viewer metadata uses these same source-header density bounds.
        output_data[mask] = to_int8(
            source.data[mask], float(source.header.dmin), float(source.header.dmax),
        )
        header = source.header.copy()

    with tempfile.TemporaryDirectory() as directory:
        output_path = Path(directory) / 'volume.mrc.gz'
        with mrcfile.new(output_path, compression='gzip') as output:
            output.set_data(output_data)
            for field in (
                'nxstart', 'nystart', 'nzstart', 'mx', 'my', 'mz',
                'cella', 'cellb', 'mapc', 'mapr', 'maps', 'ispg', 'origin',
            ):
                output.header[field] = header[field]
        compressed_data = output_path.read_bytes()

    response = HttpResponse(compressed_data, content_type='application/octet-stream')
    response['Content-Encoding'] = 'gzip'
    response['X-Content-Type-Options'] = 'nosniff'
    cache_version_matches = (
        source_version is not None
        and requested_version == source_version
        and volume_cache_version(input_file) == source_version
    )
    if cache_version_matches:
        patch_cache_control(response, private=True, max_age=VOLUME_BROWSER_CACHE_SECONDS,
                            no_transform=True)
    else:
        patch_cache_control(response, private=True, no_store=True, no_transform=True)
    return response
