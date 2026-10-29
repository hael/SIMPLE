"""Optional, bounded storage probes; never flush global caches or touch user files."""
import hashlib
import os
from pathlib import Path
import random
import shutil
import socket
import statistics
import tempfile
import time

MIB = 1024**2


def positive(value):
    import argparse
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return number


def add_options(parser):
    parser.add_argument('--io-dirs', nargs='+', type=Path,
                        help='optional storage comparison directories (must already exist)')
    parser.add_argument('--io-size-mib', type=positive, default=256,
                        help='sequential probe size per directory (default: 256 MiB)')
    parser.add_argument('--io-repeats', type=positive, default=3,
                        help='storage repetitions (default: 3)')
    parser.add_argument('--io-small-files', type=positive, default=128,
                        help='number of 64 KiB files per repetition (default: 128)')


def sync_directory(directory):
    fd = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def write_files(paths, blocks):
    started = time.perf_counter()
    for path, count in paths:
        with path.open('wb') as stream:
            for index in range(count):
                stream.write(blocks[index % len(blocks)])
            stream.flush()
            os.fsync(stream.fileno())
    sync_directory(paths[0][0].parent)
    return time.perf_counter() - started


def read_files(paths, expected):
    started = time.perf_counter()
    digests = []
    for path, _ in paths:
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            while True:
                block = stream.read(MIB)
                if not block:
                    break
                digest.update(block)
        digests.append(digest.digest())
    seconds = time.perf_counter() - started
    if digests != expected:
        raise ValueError('storage probe content verification failed')
    return seconds


def evict_hint(paths):
    """Per-file hint only: this cannot guarantee physically cold reads."""
    if not hasattr(os, 'posix_fadvise') or not hasattr(os, 'POSIX_FADV_DONTNEED'):
        return False
    try:
        for path, _ in paths:
            with path.open('rb') as stream:
                os.posix_fadvise(stream.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
        return True
    except OSError:
        return False


def benchmark(directory, size_mib=256, repeats=3, small_files=128):
    if min(size_mib, repeats, small_files) < 1:
        raise ValueError('storage sizes and repetitions must be positive')
    directory = Path(directory).resolve(strict=True)
    if not directory.is_dir():
        raise ValueError('storage target must be an existing directory')
    free = shutil.disk_usage(directory).free
    required = size_mib*MIB + small_files*65536
    if free < required + max(64*MIB, required//10):
        raise ValueError('insufficient free space for storage probe and reserve')
    # Incompressible deterministic blocks, prepared outside the timed region.
    rng = random.Random(20261029)
    blocks = [rng.randbytes(MIB) for _ in range(min(size_mib, 16))]
    small = [rng.randbytes(65536)]
    large_hash = hashlib.sha256()
    for index in range(size_mib):
        large_hash.update(blocks[index % len(blocks)])
    result = {'status': 'measured', 'directory': str(directory), 'host': socket.gethostname(),
              'device_id': directory.stat().st_dev, 'free_bytes_before': free,
              'size_mib': size_mib, 'small_files': small_files, 'small_file_bytes': 65536,
              'repetitions': repeats, 'measurements': [],
              'payload': 'deterministic pseudorandom; cycles through at most 16 MiB; deduplication may affect results',
              'scope': 'storage microbenchmark on this host; not SIMPLE workflow performance',
              'cache_policy': 'per-file eviction hint before first read; cold storage is not guaranteed',
              'write_policy': 'file fsync plus directory fsync included',
              'read_policy': 'SHA-256 verification included in read timings',
              'shared_across_nodes': None}
    with tempfile.TemporaryDirectory(prefix='.simple-io-', dir=directory) as temporary:
        root = Path(temporary)
        large_paths = [(root/'stack.bin', size_mib)]
        small_paths = [(root/('part%05d.bin' % i), 1) for i in range(small_files)]
        for repeat in range(repeats):
            for kind, paths, content, hashes, length in (
                ('sequential', large_paths, blocks, [large_hash.digest()], size_mib*MIB),
                ('small_files', small_paths, small, [hashlib.sha256(small[0]).digest()]*small_files, small_files*65536)):
                seconds = write_files(paths, content)
                rows = [('write_sync', seconds, None)]
                hinted = evict_hint(paths)
                rows.append(('first_read', read_files(paths, hashes), hinted))
                rows.append(('repeat_read', read_files(paths, hashes), False))
                for operation, elapsed, hint in rows:
                    result['measurements'].append({'repeat': repeat, 'pattern': kind, 'operation': operation,
                        'seconds': elapsed, 'mib_per_second': length/MIB/elapsed, 'bytes': length,
                        'eviction_hint_accepted': hint})
                for path, _ in paths:
                    path.unlink()
    result['summary'] = []
    for kind in ('sequential', 'small_files'):
        for operation in ('write_sync', 'first_read', 'repeat_read'):
            rows = [r for r in result['measurements'] if r['pattern'] == kind and r['operation'] == operation]
            durations = [r['seconds'] for r in rows]
            result['summary'].append({'pattern': kind, 'operation': operation,
                'median_seconds': statistics.median(durations),
                'median_mib_per_second': statistics.median(r['mib_per_second'] for r in rows),
                'timing_cv': statistics.stdev(durations)/statistics.mean(durations) if repeats > 1 else None})
    return result


def run(args):
    results = []
    for directory in args.io_dirs or []:
        print('Measuring storage: ' + str(directory), flush=True)
        try:
            result = benchmark(directory, args.io_size_mib, args.io_repeats, args.io_small_files)
        except Exception as error:
            result = {'directory': str(directory), 'host': socket.gethostname(), 'status': 'failed', 'reason': str(error)}
        results.append(result)
    return results


def report(results):
    if not results:
        return
    print('\nStorage comparison (MiB/s; writes include fsync; reads include verification)')
    print('Directory  Pattern  Operation  Median seconds  MiB/s  CV')
    for result in results:
        if result['status'] != 'measured':
            print(result['directory'] + ': ' + result['reason'])
        for row in result.get('summary', []):
            cv = '-' if row['timing_cv'] is None else f"{row['timing_cv']:.1%}"
            print(f"{result['directory']}  {row['pattern']}  {row['operation']}  {row['median_seconds']:.4f}  {row['median_mib_per_second']:.1f}  {cv}")
    print('First reads use an eviction hint, not a guaranteed cold cache. Node-local/shared storage is not inferred.')
