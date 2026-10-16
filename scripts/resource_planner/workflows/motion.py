"""Motion fixture, truth validation and scientific settings."""
from array import array
import math
import statistics
import struct
import sys
from memory.benchmark_motion_correct import make_mrc_header, synthetic_row

SHIFTS = [0, 0, 0, 0, 0, 1, 2, 3]
SCIENCE = {'algorithm':'patch','mcpatch':'yes','downscale':'yes','smpd_downscale':1.3}

def generate(root, size):
    root.mkdir()
    base, movie = root/'truth.mrc', root/'movie.mrc'
    with base.open('wb') as f:
        f.write(make_mrc_header(size, size, 1, 1.3))
        for y in range(size):
            f.write(synthetic_row(size, size, y))
    with movie.open('wb') as out, base.open('rb') as source:
        out.write(make_mrc_header(size, size, 8, 1.3))
        for shift in SHIFTS:
            source.seek(1024)
            for y in range(size):
                row = source.read(size*4)
                if shift:
                    row = row[-shift*4:] + row[:-shift*4]
                out.write(row)
    return movie


def image(path, size):
    with path.open('rb') as f:
        header = f.read(1024)
        if struct.unpack_from('<4i', header) != (size, size, 1, 2):
            raise ValueError(f"unexpected image header: {path}")
        f.seek(1024 + struct.unpack_from('<i', header, 92)[0])
        data = array('f')
        data.fromfile(f, size*size)
    if sys.byteorder != 'little':
        data.byteswap()
    if not all(math.isfinite(v) for v in data):
        raise ValueError('nonfinite image')
    return data


def validate(directory, truth, size):
    data = image(directory/'movie_intg.mrc', size)
    step = max(1, size*size//65536)
    x, y = truth[::step], data[::step]
    mx, my = statistics.mean(x), statistics.mean(y)
    correlation = sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(
        sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y))
    rows = []
    for line in (directory/'movie.star').read_text().splitlines():
        fields = line.split()
        if len(fields) == 4 and fields[0].isdigit():
            rows.append([float(v) for v in fields])
    if len(rows) != 8 or not all(math.isfinite(v) for row in rows for v in row):
        raise ValueError('invalid global shift table')
    # Image fidelity is the gate; retain the full shift estimate for review.
    return {'passed': correlation >= 0.98, 'correlation': correlation,
            'shifts_x': [r[1] for r in rows], 'shifts_y': [r[2] for r in rows]}



from resource_planner.workflows.common import ImageAdapter
class MotionAdapter(ImageAdapter):
    name = 'motion_correct'
    project_name = 'motion'
    rate = 'movies_per_second'
    science = SCIENCE
    workload = {'frames': 8, 'smpd': 1.3, 'algorithm': 'patch', 'mcpatch': 'yes'}
    validator = staticmethod(validate)

    def metadata(self, args):
        return {**super().metadata(args), **self.workload, 'known_shifts_x': SHIFTS}


