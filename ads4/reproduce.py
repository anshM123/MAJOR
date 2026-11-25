#!/usr/bin/env python3
"""reproduce.py: one command for the checks behind the headline claims of the Letter.

    python3 reproduce.py           # quick mode (a few CPU-minutes)
    python3 reproduce.py --full    # adds the T0 certificate (~20 CPU-min) and the longer proof scripts
    python3 reproduce.py --list    # list the stages without running them

Each stage runs an existing project script (unchanged), checks its exit status and, where given, that its
output contains an expected string; a summary table is printed at the end and written to reproduce_report.txt.
Requirements: python3 with numpy, scipy, sympy, mpmath, matplotlib, python-flint (Arb ball arithmetic).
Every script writes only into its own cache/ directory (or prl/fig/ for the figures).
"""
import os, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))
ENV = dict(os.environ, OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')

# (name, claim, working dir, command, expected substrings, full-only)
STAGES = [
    ("tau0", "Result (1): tau_0 at K=1 from the published weights; Fig. 1",
     "prl/fig", ["python3", "make_fig.py"], ["tau0 = 2.89126e-08"], False),
    ("fig2", "Fig. 2 (coverage map)", "prl/fig", ["python3", "make_fig2.py"], [], False),
    ("K", "Result (1): certified K = 1 - 32.4/Lambda^2 inputs (delta K0, grid, bounds, self-test)",
     "code/kcert", ["bash", "kc_run_all.sh"], [], False),
    ("H17", "Result (1): closed-form H17, certified K and tau_0 <= 2.8914e-8",
     "code/kb2a4", ["python3", "-B", "kb_h17.py"], [], False),
    ("ang-cdh", "Angular route: celestial recurrence = continuous dual Hahn; K_B2A4 polynomial",
     "code/angular", ["python3", "verify_cdh.py"], ["OK"], False),
    ("ang-ward", "Angular route: all-J first-moment operator (rest-frame Ward algebra)",
     "code/angular", ["python3", "verify_ward.py"], ["diagonal T+C-2m^2 OK"], False),
    ("ang-ladder", "Angular route: all-J helicity ladder (closed forms, normalised ratio)",
     "code/angular", ["python3", "verify_ladder.py"], ["OK"], False),
    ("ang-sp3", "Angular route: SP3 positivity of the first-moment operator; smeared coefficients",
     "code/angular", ["python3", "verify_sp3.py"], ["OK"], False),
    ("ang-s6", "Angular route: generating function, Y0 = -Einstein, Y4 = 303/5120, gap constant",
     "code/angular", ["python3", "verify_s6.py"], ["Y_4 = 303/5120"], False),
    ("ang-s9", "Angular route: telescoping and exact endpoint identity",
     "code/angular", ["python3", "verify_s9.py"], [], False),
    ("ang-s10", "Angular route: coercivity, contact kernel Q(theta), K(g) and tau <= 0.0178",
     "code/angular", ["python3", "verify_s10.py"], ["0.01781"], False),
    ("ang-s13", "Angular route: truncation identities of the defect theorem",
     "code/angular", ["python3", "verify_s13.py"], ["1/24"], False),
    ("ang-const", "Angular route: uniform constant K(g), fourth-harmonic ratio (independent)",
     "code/angmatch", ["python3", "am2_constant.py"], [], False),
    ("ang-ward2", "Angular route: Ward step C(0) = C(1) (independent)",
     "code/angmatch", ["python3", "am3_ward.py"], [], False),
    # ---------------------------------------------------------------- full mode
    ("T0", "T0: certified leading-Regge floor >= 0.1999 N_s (Delta_g >= 1000), 0.199 (>= 300)",
     "code/t0cert", ["bash", "run_all.sh"], ["VERDICT"], True),
    ("B2A4", "K_B2A4 in full (one-point traces, k=2 ratio, k=3 celestial polynomial, H17)",
     "code/kb2a4", ["bash", "kb_run_all.sh"], [], True),
    ("collinear", "Crossing proof: all-J conserved structures and collinear identities",
     "code/angproofs", ["python3", "ap3_collinear.py"], [], True),
    ("extreme", "Crossing proof: extreme detector amplitudes (Chu-Vandermonde)",
     "code/angproofs", ["python3", "ap4_extreme.py"], [], True),
    ("envelopes", "Crossing proof: kernel envelopes with explicit constants (Wendel)",
     "code/angproofs", ["python3", "ap6_envelopes.py"], [], True),
    ("blocks", "Crossing proof: collinear block lower bound (Euler integral)",
     "code/angproofs", ["python3", "ap8_blocks.py"], [], True),
    ("lowspin", "Crossing proof: low-spin sectors",
     "code/angproofs", ["python3", "ap7_lowspin.py"], [], True),
    ("free", "Free CFT3s: null-plane algebra and endpoint rule 3/2 + 3/2 = 3",
     "code/freenp", ["python3", "np5_spectral.py"], [], True),
]


def main():
    full = '--full' in sys.argv
    if '--list' in sys.argv:
        for name, claim, d, cmd, exp, fo in STAGES:
            print(f"{'[full] ' if fo else '       '}{name:10s} {claim}\n           {d}: {' '.join(cmd)}")
        return 0
    rows, t_all = [], time.time()
    for name, claim, d, cmd, exp, fo in STAGES:
        if fo and not full:
            continue
        wd = os.path.join(ROOT, d)
        os.makedirs(os.path.join(wd, 'cache'), exist_ok=True)
        t = time.time()
        print(f"== {name}: {claim}\n   ({d}: {' '.join(cmd)})", flush=True)
        try:
            p = subprocess.run(cmd, cwd=wd, env=ENV, capture_output=True, text=True, timeout=7200)
            out = p.stdout + p.stderr
            ok = p.returncode == 0 and all(e in out for e in exp)
            why = '' if ok else (f'exit {p.returncode}' if p.returncode else 'expected output missing: ' + ', '.join(e for e in exp if e not in out))
        except subprocess.TimeoutExpired:
            ok, out, why = False, '', 'timeout'
        with open(os.path.join(wd, 'cache', f'reproduce_{name}.log'), 'w') as f:
            f.write(out)
        dt = time.time() - t
        rows.append((name, 'PASS' if ok else 'FAIL', dt, why, claim))
        print(f"   -> {'PASS' if ok else 'FAIL ' + why}  [{dt:.0f}s]", flush=True)
    lines = ['', f"reproduce.py ({'full' if full else 'quick'} mode), {time.strftime('%Y-%m-%d %H:%M')} UTC, "
             f"total {time.time() - t_all:.0f}s", '']
    lines += [f"{st:4s}  {name:10s} {dt:6.0f}s  {claim}" + (f"   [{why}]" if why else '') for name, st, dt, why, claim in rows]
    nfail = sum(r[1] != 'PASS' for r in rows)
    lines += ['', f"{len(rows) - nfail}/{len(rows)} stages passed; logs in <dir>/cache/reproduce_<stage>.log"]
    print('\n'.join(lines))
    with open(os.path.join(ROOT, 'reproduce_report.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
