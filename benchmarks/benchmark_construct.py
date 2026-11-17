"""Run isolated native mesh construction benchmarks; no Python package install needed."""

import os
import csv
import sys
import json
import time
import hashlib
import argparse
import platform
import statistics
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_optional(path):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def git_output(*args):
    try:
        return subprocess.check_output(['git', '-C', str(ROOT), *args], text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def metadata(binary, timeout):
    result = subprocess.run([str(binary), '--info'], capture_output=True, text=True, check=True, timeout=timeout)
    cache = {}
    for parent in binary.parents:
        path = parent / 'CMakeCache.txt'
        if path.is_file():
            for line in path.read_text().splitlines():
                if '=' in line and not line.startswith(('#', '//')):
                    key, value = line.split('=', 1)
                    key = key.split(':', 1)[0]
                    if key.startswith(('CMAKE_CXX_', 'CMAKE_BUILD_TYPE', 'TINYMESH_')):
                        cache[key] = value
            break
    cpu = platform.processor()
    cpuinfo = read_optional('/proc/cpuinfo')
    if cpuinfo:
        cpu = next(
            (line.split(':', 1)[1].strip() for line in cpuinfo.splitlines() if line.startswith('model name')), cpu
        )
    status = git_output('status', '--porcelain')
    return {
        'timestamp_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'platform': platform.platform(),
        'cpu': cpu,
        'logical_cpus': os.cpu_count(),
        'affinity_cpus': len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
        'cgroup_cpu_max': read_optional('/sys/fs/cgroup/cpu.max'),
        'cgroup_memory_max': read_optional('/sys/fs/cgroup/memory.max'),
        'python': sys.version,
        'binary': str(binary),
        'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
        'source_revision': git_output('rev-parse', 'HEAD'),
        'source_dirty': status is not None and bool(status),
        'mesh_cpp_sha256': hashlib.sha256((ROOT / 'src/tinymesh/core/mesh.cpp').read_bytes()).hexdigest(),
        'native_build': json.loads(result.stdout),
        'cmake_cache': cache,
        'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS'),
    }


def measure(binary, job, seed, timeout):
    if job['file']:
        command = [str(binary), '--file', job['file']]
    else:
        command = [
            str(binary),
            '--case',
            job['case'],
            '--size',
            str(job['size']),
            '--order',
            job['order'],
            '--seed',
            str(seed),
        ]
    start = time.perf_counter()
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {'status': 'timeout', 'wall_seconds': time.perf_counter() - start, 'error': f'Exceeded {timeout}s'}
    row = {'wall_seconds': time.perf_counter() - start, 'returncode': result.returncode}
    if result.returncode:
        return {**row, 'status': 'error', 'error': (result.stderr + result.stdout)[-4000:]}
    try:
        # Mesh readers may print warnings before the final JSON record.
        record = json.loads(result.stdout.strip().splitlines()[-1])
        if not record['valid'] or record['build_seconds'] <= 0:
            raise ValueError('Invalid construction result or timer')
        return {**row, **record, 'status': 'ok', 'stderr': result.stderr}
    except (ValueError, KeyError, IndexError) as error:
        return {**row, 'status': 'error', 'error': f'{error}: {result.stdout[-4000:]}'}


def summarize(samples):
    groups = {}
    for row in samples:
        if not row['warmup']:
            key = tuple(row[k] for k in ('mode', 'case', 'size', 'order', 'file'))
            groups.setdefault(key, []).append(row)
    summaries = []
    for key, rows in groups.items():
        good = [r for r in rows if r['status'] == 'ok']
        summary = dict(zip(('mode', 'case', 'size', 'order', 'file'), key, strict=True))
        summary.update(successful_runs=len(good), failed_runs=len(rows) - len(good))
        if good:
            values = [r['build_seconds'] for r in good]
            summary.update(
                triangles=good[0]['triangles'],
                vertices=good[0]['vertices'],
                median_seconds=statistics.median(values),
                min_seconds=min(values),
                max_seconds=max(values),
                stdev_seconds=statistics.stdev(values) if len(values) > 1 else 0.0,
                median_triangles_per_second=statistics.median(r['triangles'] / r['build_seconds'] for r in good),
                median_preparation_seconds=statistics.median(r['preparation_seconds'] for r in good),
                median_verification_seconds=statistics.median(r['verification_seconds'] for r in good),
                median_destruction_seconds=statistics.median(r['destruction_seconds'] for r in good),
                max_peak_rss_bytes=max(
                    (r['peak_rss_bytes'] for r in good if r['peak_rss_bytes'] is not None), default=None
                ),
            )
        summaries.append(summary)
    return summaries


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--cases', nargs='+', choices=['grid', 'torus', 'triangles', 'fan'])
    parser.add_argument('--sizes', nargs='+', type=int, default=[128, 256, 512])
    parser.add_argument('--orders', nargs='+', choices=['ordered', 'shuffled'], default=['ordered'])
    parser.add_argument('--files', nargs='+', type=Path, default=[])
    parser.add_argument('--repeat', type=int, default=3)
    parser.add_argument('--warmup', type=int, default=1)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--timeout', type=float, default=120)
    parser.add_argument('--output', type=Path, default=ROOT / 'benchmarks/results/construct')
    args = parser.parse_args()
    if args.repeat < 1 or args.warmup < 0 or args.timeout <= 0 or min(args.sizes) < 3:
        parser.error('repeat >= 1, warmup >= 0, timeout > 0, and sizes >= 3 are required')
    if not 0 <= args.seed <= 2**32 - 1:
        parser.error('seed must be an unsigned 32-bit integer')
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f'Benchmark executable not found: {binary}')

    cases = args.cases if args.cases is not None else ([] if args.files else ['grid', 'torus'])
    jobs = [
        {'mode': 'construct', 'case': case, 'size': size, 'order': order, 'file': ''}
        for case in cases
        for size in args.sizes
        for order in args.orders
    ]
    jobs += [{'mode': 'load', 'case': '', 'size': None, 'order': '', 'file': str(p.resolve())} for p in args.files]
    environment = metadata(binary, args.timeout)
    if environment['native_build']['build_config'] not in {'Release', 'RelWithDebInfo'}:
        print('Warning: use a Release or RelWithDebInfo executable for performance comparisons.', file=sys.stderr)
    samples = []
    prefix = args.output.resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    failed = False
    for job in jobs:
        label = job['file'] or f'{job["case"]}/{job["size"]}/{job["order"]}'
        for iteration in range(-args.warmup, args.repeat):
            sample = {**job, 'warmup': iteration < 0, 'iteration': iteration, 'seed': args.seed}
            sample.update(measure(binary, job, args.seed, args.timeout))
            samples.append(sample)
            failed |= sample['status'] != 'ok'
            timing = f'{sample["build_seconds"]:.6f}s' if sample['status'] == 'ok' else sample['status']
            print(f'{label} run={iteration}: {timing}', file=sys.stderr, flush=True)
        # Keep completed samples if a later case is interrupted.
        report = {
            'metadata': environment,
            'settings': vars(args) | {'binary': str(binary), 'output': str(prefix)},
            'samples': samples,
            'summary': summarize(samples),
        }
        report['settings']['files'] = [str(p.resolve()) for p in args.files]
        Path(f'{prefix}.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
        write_csv(Path(f'{prefix}.samples.csv'), samples)
        write_csv(Path(f'{prefix}.summary.csv'), report['summary'])
    print(f'Results: {prefix}.json', flush=True)
    return int(failed)


if __name__ == '__main__':
    raise SystemExit(main())
