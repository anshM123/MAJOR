"""truncate_grids.py -- truncated copies of the stored k = 2 Abel grid weights (v8 of the Supplement; NU_TAIL_CHECK.md, problem (c)).

Why.  The stored minimal-norm grid weights d_abel_X_{300,1000}.txt (Omega^EE) and k2/xq_abel_X_{300,1000}.txt (X^Q of Omega^Q_up)
extend past Lambda = sqrt((Dg-2)(Dg-1)): to Lambda + 29 (Dg = 300) and Lambda + 312 (Dg = 1000).  Beyond Lambda + 10 they contain
no information about the weight: the exact weight has decayed like exp(-(pi/2)(nu - Lambda)) to below 1e-15 there, and what is left
in the files is a double-precision floor (|X| ~ 1e-16 .. 1e-14) that does not decay.  Read literally against the exact per-nu image
of a near-gap heavy state, which grows exponentially once nu exceeds the state's mass, that floor would contribute an uncontrolled
tail at Dg = 1000 (NU_TAIL_CHECK Sec. 3.2).  The Supplement (Sec. S2.B-C) therefore DEFINES the k = 2 pieces by the continuum
weights k2/k2_continuum_{D,Q}_*.json and states that the grids tabulate them on nu <= Lambda + 10 and are to be read as zero beyond.
No computed number uses the grids beyond Lambda (every per-state evaluation is on nu <= Lambda <= m), so this changes no number.

What.  For each of the four grids, write a copy that keeps the header and the rows with nu_k <= Lambda + 10 (the rows beyond are
dropped, i.e. read as zero).  The published files are NOT modified (their hashes are quoted in MANIFEST.txt and k2/SHA256SUMS).
Outputs (new files): d_abel_X_300_trunc.txt, d_abel_X_1000_trunc.txt, k2/xq_abel_X_300_trunc.txt, k2/xq_abel_X_1000_trunc.txt,
and TRUNC_SHA256SUMS (full SHA-256 of the four copies).  Printed checks per grid: number of rows kept / dropped, the largest |X|
dropped, the dropped L1 mass sum |X| dnu (absolute and relative to the kept mass), and, as a consistency check, the largest
|X_grid - X_continuum| on (Lambda, Lambda + 10] from the continuum tables (same nodes).

Standalone (Python 3 standard library only).  Usage (from ads4/prl/data):  python3 truncate_grids.py
"""
import hashlib
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MARGIN = 10.0                      # rows with nu_k > Lambda + MARGIN are dropped
GRIDS = [  # (source, continuum table for the check, Delta_g, dnu)
    ('d_abel_X_300.txt', 'k2/k2_continuum_D_300_table.txt', 300, 0.02),
    ('d_abel_X_1000.txt', 'k2/k2_continuum_D_1000_table.txt', 1000, 0.04),
    ('k2/xq_abel_X_300.txt', 'k2/k2_continuum_Q_300_table.txt', 300, 0.02),
    ('k2/xq_abel_X_1000.txt', 'k2/k2_continuum_Q_1000_table.txt', 1000, 0.04),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_x(path):
    """{nu string: (Re X, Im X)} from a table whose columns 4, 5 are Re X, Im X (continuum tables)."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.split()
            out[c[0]] = (float(c[3]), float(c[4]))
    return out


def main():
    sums = []
    for src, cont, dg, dnu in GRIDS:
        lam = math.sqrt((dg - 2) * (dg - 1))
        cut = lam + MARGIN
        p_src = os.path.join(HERE, src)
        p_out = p_src[:-4] + '_trunc.txt'
        xc = read_x(os.path.join(HERE, cont))
        header, kept, n_drop, started = [], 0, 0, False
        drop_max, drop_l1, kept_l1, cmp_max, psi_beyond = 0.0, 0.0, 0.0, 0.0, 0.0
        with open(p_src) as f, open(p_out, 'w') as g:
            for line in f:
                if line.startswith('#'):
                    header.append(line)
                    continue
                if not started:                      # first data row: write the header with the provenance line
                    started = True
                    g.write(header[0])
                    g.write('# TRUNCATED COPY (truncate_grids.py, Supplement v8): rows with nu_k > Lam + %g = %.4f removed, i.e. the '
                            'weight is read as zero there; source %s, sha256 %s\n' % (MARGIN, cut, os.path.basename(src),
                                                                                     sha256(p_src)))
                    for h in header[1:]:
                        g.write(h)
                c = line.split()
                nu = float(c[0]); xr, xi = float(c[1]), float(c[2]); ax = math.hypot(xr, xi)
                if nu > lam:
                    psi_beyond = max(psi_beyond, abs(float(c[5])))
                if nu <= cut:
                    g.write(line); kept += 1; kept_l1 += ax * dnu
                    if lam < nu and c[0] in xc:
                        cr, ci = xc[c[0]]
                        cmp_max = max(cmp_max, math.hypot(xr - cr, xi - ci))
                else:
                    n_drop += 1; drop_max = max(drop_max, ax); drop_l1 += ax * dnu
        sums.append((sha256(p_out), os.path.relpath(p_out, HERE)))
        print('%-24s Dg=%4d Lam=%.4f  kept %6d rows (nu <= %.2f), dropped %6d;  max|X| dropped %.2e;  sum|X|dnu dropped %.2e '
              '(%.1e of kept);  max|X_grid - X_cont| on (Lam, Lam+10] %.1e;  max|psi2| beyond Lam %.1e'
              % (src, dg, lam, kept, cut, n_drop, drop_max, drop_l1, drop_l1 / kept_l1, cmp_max, psi_beyond))
    with open(os.path.join(HERE, 'TRUNC_SHA256SUMS'), 'w') as g:
        for h, p in sums:
            g.write('%s  %s\n' % (h, p))
    print('wrote the four *_trunc.txt copies and TRUNC_SHA256SUMS; published files unchanged')


if __name__ == '__main__':
    main()
