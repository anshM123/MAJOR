"""
ext/run_ext.py -- run every job in jobs.txt on N worker processes (default: all cores minus 2).
Jobs already finished (their 'done' marker exists) are skipped, so the run can be stopped and restarted at any time.
Usage: python3 run_ext.py [N] [jobs.txt]
Then:  python3 pack.py   -> results/ext_results.json  (commit and push that one file)
"""
import os, sys, subprocess, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
N = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, (os.cpu_count() or 4) - 2)
JOBS = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'jobs.txt')


def parse():
    jobs = []
    for line in open(JOBS):
        line = line.split('#')[0].split()
        if not line:
            continue
        tag, meth, J, D, seed = line[0], line[1], int(line[2]), float(line[3]), int(line[4])
        for nu in line[5].split(','):
            jobs.append((tag, meth, J, D, float(nu), seed))
    return jobs


def one(job):
    tag, meth, J, D, nu, seed = job
    out = os.path.join(HERE, 'out', tag, f'{meth}_J{J}_D{D:g}_nu{nu:.6f}_s{seed}')
    done = os.path.join(out, 'done')
    if os.path.exists(done):
        return job, 'skip', 0.0
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(out, 'log.txt'), 'w') as lg:
        rc = subprocess.call([sys.executable, os.path.join(HERE, 'point.py'), meth, str(J), repr(D), repr(nu), str(seed), out],
                             stdout=lg, stderr=subprocess.STDOUT)
    if rc == 0:
        open(done, 'w').write('ok\n')
    return job, ('ok' if rc == 0 else f'FAILED rc={rc}'), time.time() - t0


if __name__ == '__main__':
    jobs = parse()
    print(f'{len(jobs)} jobs on {N} workers', flush=True)
    with Pool(N) as p:
        for i, (job, st, dt) in enumerate(p.imap_unordered(one, jobs), 1):
            print(f'[{i}/{len(jobs)}] {st:8s} {dt/60:6.1f} min  {job}', flush=True)
