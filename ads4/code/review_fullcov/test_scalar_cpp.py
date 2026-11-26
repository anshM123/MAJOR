"""PROOF.md §4/§6 code faithfulness: run the portable GMP scalar (scalar_ref_portable.cpp,
compiled into the scratch dir) on a few nodes taken from archived shard inputs, and compare
with the independent fixed-point evaluator rc.F_K at the same rational (c', h').
Checks: |scalar F - exact F_K(c',h')| <= E_scalar, the node-y convention y'=h'/H(c'),
and that the printed [lo,hi] encloses F_K(c',h') - R(x',y') up to the reference error."""
import sys, subprocess, os
from fractions import Fraction as Fr
from pathlib import Path
import rc

INDIR = Path(os.environ.get('SHARD_DIR', '/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/fullcov/sliver_merge'))
SCR = Path(os.environ.get('REVIEW_SCRATCH', '/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/review_fullcov'))
SCR.mkdir(parents=True, exist_ok=True)
exe = SCR / 'scalar_ref_portable'
if not exe.exists():
    subprocess.run(['g++', '-O2', '-std=c++17', str(rc.MERGE / 'scalar_ref_portable.cpp'), '-lgmpxx', '-lgmp', '-o', str(exe)], check=True)

def mini(shard, rows=1, ys=(0, 1, -1)):
    tok = (INDIR / shard).read_text().split('\n')
    K, Den, Ny, Nc, P = map(int, tok[0].split())
    shifts = tok[1]; yden = tok[2]; Y = tok[3].split(); refs = tok[4]; R = tok[5:5 + Nc]
    Ysel = [Y[i] for i in ys]
    txt = f'{K} {Den} {len(Ysel)} {rows} {P}\n{shifts}\n{yden}\n{" ".join(Ysel)}\n{refs}\n' + '\n'.join(R[:rows]) + '\n'
    return txt, K, Den, P, int(yden), list(map(int, Ysel)), list(map(int, refs.split())), [list(map(int, r.split())) for r in R[:rows]]

def check(stem, shard, b_list=(0, 4)):
    z, l, u, Hl, Hu = rc.cell(stem); Kfull = z['K']
    txt, K, Den, P, Yden, Ysel, refs, rows = mini(shard)
    res = subprocess.run([str(exe), rc.COEFF, 'debug'], input=txt, capture_output=True, text=True, check=True)
    dbg = [ln.split() for ln in res.stderr.splitlines() if len(ln.split()) == 4]
    out = {int(x.split()[0]): tuple(map(int, x.split()[1:])) for x in res.stdout.splitlines() if not x.startswith('COUNT')}
    err = (Kfull + 1) * (Fr(4) + Hl * ((Kfull + 1) ** 2 + Kfull * Kfull + 8 * Kfull + 67)) / 10 ** 60 + Fr(10 * (Kfull + 1), 2 ** P)
    ln, ld, un, ud = refs[:4]; ref = refs[4:]
    worst = 0; enc = True; mins = {}; maxs = {}
    for ic, j, b, v in dbg:
        ic, j, b, v = int(ic), int(j), int(b), int(v)
        C, Hn, Hd = rows[ic]; c = Fr(C, Den); Hc = Fr(Hn, Hd)
        assert Hc == rc.Hsec(c, l, u, Hl, Hu)
        hint = (Hn * Den * Ysel[j]) // (Hd * Yden); h = Fr(hint, Den)
        Fsc = Fr(v, 2 ** P) * (1 - c) ** 3
        Fex = rc.F_K(c, h, K, b, W=300)
        d = abs(Fsc - Fex); worst = max(worst, d / err)
        x = (c - l) / (u - l); y = h / Hc
        Pm = [ref[(0 if b == 0 else 9) + 3 * a_ + d_] for a_ in range(3) for d_ in range(3)]
        Pm = [Pm[0:3], Pm[3:6], Pm[6:9]]
        g = Fex - rc.Rval(Pm, x, y)
        mins[b] = min(mins.get(b, g), g); maxs[b] = max(maxs.get(b, g), g)
    referr = Fr(4, 1000 * 2 ** 50) + Fr(8, 10 ** 20)
    for b in b_list:
        lo, hi = out[b]
        ok = Fr(lo, 10 ** 20) - err - referr <= mins[b] and maxs[b] <= Fr(hi, 10 ** 20) + err + referr
        enc &= ok
        print(f'  {shard} b={b}: printed [lo,hi]=[{lo},{hi}]e-20  independent F_K-R in [{float(mins[b]):.6e},{float(maxs[b]):.6e}]  {"encloses" if ok else "MISMATCH"}')
    print(f'{shard}: K={K} P={P} nodes={len(dbg)}  max |scalar F - F_K| / E_scalar = {float(worst):.3e}', 'OK' if worst <= 1 and enc else 'FAIL', flush=True)

if __name__ == '__main__':
    check('cell_8936_9000', 'cell_8936_9000_r0.in')
    check('cell_5305_6041', 'cell_5305_6041_g3_0.in')
    check('cell_7596_7929', 'cell_7596_7929_s3_1.in')
