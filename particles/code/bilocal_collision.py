#!/usr/bin/env python3
"""The explicit two-massless-particle interval kernel and its matching model.

Source: Aref'eva, Khramtsov, Tikhanovskaya, arXiv:1706.07390, eqs. (3.40),
(4.5), and (4.6). Length units are 4 G_N S = 6 S/c, with AdS radius one.
This program DOES NOT establish simultaneous HRT admissibility for unions.
Its multiparty output is the noncrossing-matching functional of the published
single-interval kernel. See COLLISION_KERNEL.md for the precise theorem.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import json
from pathlib import Path
import math
import numpy as np

PI = math.pi


def interval_length(a: float, b: float, t: float, radius: float = 1.0,
                    epsilon: float = 1e-4) -> float:
    """Symmetric pure-state interval kernel on the circle of circumference 2pi.

    The BTZ rest-frame time t>=0 is NOT global collision time. At a finite
    t the infalling particles are at r=R coth(Rt), outside the final horizon.
    epsilon is our dimensionless endpoint cutoff, used consistently in both
    branches; it is not the paper's ambiguous squared-cutoff notation.
    """
    if t < 0 or radius <= 0 or epsilon <= 0:
        raise ValueError("Require t>=0, R>0, epsilon>0")
    delta = (b - a) % (2 * PI)
    if delta < 1e-13 or 2 * PI - delta < 1e-13:
        return 0.0
    if delta > PI:
        a, delta = b, 2 * PI - delta
    start = a % (2 * PI)
    end = start + delta
    d = radius * delta / 2
    direct = math.log(4 * math.sinh(d) ** 2)
    p = (math.floor(start / PI) + 1) * PI
    # A shortest arc contains at most one injection point, 0 or pi.
    if start + 1e-13 < p < end - 1e-13:
        center = radius * ((start + end) / 2 - p)
        E = 1 / math.tanh(PI * radius / 2)
        P = math.sinh(d) - E * math.cosh(d)
        Q = E * math.cosh(radius * t) + math.cosh(center) * P
        # Algebraically equal to eq. (3.40), with less cancellation.
        argument = 4 * (Q * Q + math.sinh(center) ** 2 * (E * E - P * P))
        if argument <= 0:
            raise ArithmeticError((a, b, t, radius, argument))
        direct = min(direct, math.log(argument))
    return direct + 2 * math.log(1 / epsilon)


def kernel(endpoints, t, radius=1.0, epsilon=1e-4):
    x = np.asarray(endpoints, dtype=float)
    return np.array([[interval_length(a, b, t, radius, epsilon) if i != j else 0.0
                      for j, b in enumerate(x)] for i, a in enumerate(x)])


def split_weights(D):
    n = len(D)
    w = np.zeros_like(D)
    for i in range(n):
        for j in range(i + 1, n):
            value = (D[i, j] + D[(i+1) % n, (j+1) % n]
                     - D[i, (j+1) % n] - D[(i+1) % n, j]) / 2
            w[i, j] = w[j, i] = value
    return w


def reconstruct(w):
    n = len(w)
    D = np.zeros_like(w)
    for a in range(n):
        for b in range(a + 1, n):
            D[a, b] = D[b, a] = sum(w[i, j] for i in range(a, b)
                for j in list(range(b, n)) + list(range(a)))
    return D


def matching_entropy(D, selected_parties):
    points = tuple(sorted(k for p in selected_parties for k in (2*p, 2*p+1)))

    @lru_cache(None)
    def solve(xx):
        if not xx:
            return 0.0
        return min(D[xx[0], xx[j]] + solve(xx[1:j]) + solve(xx[j+1:])
                   for j in range(1, len(xx), 2))
    return solve(points)


def entropy_vector(D, subset_order, parties=5):
    return np.array([matching_entropy(D, tuple(ord(c)-ord('A') for c in subset
                                               if ord(c)-ord('A') < parties))
                     for subset in subset_order])


def algebra_checks():
    """Independent symbolic identities, plus numeric physical-domain checks."""
    import sympy as s
    A, B, T, E = s.symbols('A B T E', positive=True)
    F = E*(A-T)*(B-T) + T*(B-A)
    G = E*(A*T-1)*(B*T-1) + T*(B-A)
    assert s.simplify(F*s.diff(F, A, B)-s.diff(F, A)*s.diff(F, B)-T*T) == 0
    assert s.simplify(G*s.diff(G, A, B)-s.diff(G, A)*s.diff(G, B)-T*T) == 0
    for f, da, db in ((F, (B-T)**2, -(A-T)**2),
                      (G, (B*T-1)**2, -(A*T-1)**2)):
        r = f/(T*(B-A))
        assert s.simplify(s.diff(r,A)-E*da/(T*(B-A)**2)) == 0
        assert s.simplify(s.diff(r,B)-E*db/(T*(B-A)**2)) == 0
    cd=(B/A+A/B)/2
    c2c=(A*B+1/(A*B))/2
    c2t=(T*T+1/(T*T))/2
    sd=(B/A-A/B)/2
    term=(B-1/B-A+1/A)/4-E*(B+1/B+A+1/A)/4
    published=2*(-1+E*E+(1+E*E)*cd+E*E*c2c+E*E*c2t
                 +4*E*(T+1/T)/2*term-2*E*sd)
    assert s.simplify(published-F*G/(A*B*T*T)) == 0
    # Symmetric intervals switch exactly at t=half angular width.
    worst = 0.0
    for R in (.4, 1., 2.):
        for half in (.12, .4, 1., 1.5):
            observed = interval_length(-half, half, half, R, 1.0)
            expected = math.log(4*math.sinh(R*half)**2)
            worst = max(worst, abs(observed-expected))
    assert worst < 1e-11
    return {'factorization': 'exact symbolic PASS',
            'bilinear_determinants': 'exact symbolic PASS',
            'endpoint_derivatives': 'exact symbolic PASS',
            'symmetric_switch_max_error': worst}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--facets', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('collision_results'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = json.loads(args.facets.read_text())
    facets = np.asarray(data['facets'], dtype=np.int64)
    order = data['subset_order']
    assert facets.shape == (372,31)
    algebra = algebra_checks()
    rng = np.random.default_rng(20261128)
    example = np.array([-3.02,-2.51,-2.18,-1.34,-1.10,.26,.53,1.20,1.51,2.83])
    geometries = [example]
    while len(geometries) < 8:
        x = np.sort(rng.uniform(-PI, PI, 10))
        if np.min(np.diff(np.r_[x, x[0]+2*PI])) > .045:
            geometries.append(x)
    times = np.linspace(0., 2., 81)
    all_min = float('inf')
    weight_min = float('inf')
    reconstruction_error = 0.0
    orbit_min = np.full(len(data['orbit_sizes']), np.inf)
    # Rows are globally sorted, not grouped in consecutive orbits.
    orbit_rows = [np.asarray(ii, dtype=int) for ii in data['facet_orbits']]
    assert [len(ii) for ii in orbit_rows] == data['orbit_sizes']
    vectors = 0
    worst_case = None
    curves = []
    for gi, x in enumerate(geometries):
        for R in (.5, 1., 2.):
            for t in times:
                D = kernel(x, float(t), R)
                w = split_weights(D)
                weight_min = min(weight_min, float(w[np.triu_indices(10,1)].min()))
                reconstruction_error = max(reconstruction_error,
                                           float(np.max(np.abs(reconstruct(w)-D))))
                for n in (3,4,5):
                    sv = entropy_vector(D, order, n)
                    residuals = facets @ sv
                    this_min = float(residuals.min())
                    if this_min < all_min:
                        all_min = this_min
                        worst_case = {'geometry': gi, 'R': R, 'time': float(t),
                                      'parties': n, 'facet': int(residuals.argmin())}
                    for oi in range(len(orbit_min)):
                        orbit_min[oi] = min(orbit_min[oi],
                            float(residuals[orbit_rows[oi]].min()))
                    vectors += 1
                    if gi==0 and R==1. and n==5:
                        # UV-finite I3 and cyclic inequality.
                        lookup = dict(zip(order,sv))
                        i3 = lookup['A']+lookup['B']+lookup['C']-lookup['AB']-lookup['AC']-lookup['BC']+lookup['ABC']
                        cyc = sum(lookup[''.join(sorted('ABCDE'[(i+j)%5] for j in range(3)))]
                                  -lookup[''.join(sorted('ABCDE'[(i+j)%5] for j in range(2)))] for i in range(5))-lookup['ABCDE']
                        curves.append({'t':float(t), 'length_C_ren':float(D[4,5]-2*math.log(1e4)),
                                       'minus_I3':float(-i3),'cyclic5':float(cyc),
                                       'minimum_all372':this_min,
                                       'entropies_length_units':sv.tolist()})
    assert weight_min > -1e-10
    assert all_min > -1e-10
    assert reconstruction_error < 1e-10
    report = {
        'status':'PASS for the matching model; physical union admissibility remains unproved',
        'source':'https://arxiv.org/abs/1706.07390',
        'units':'length = 4 G_N S = 6 S/c; AdS radius one',
        'time':'BTZ rest-frame boundary time, not global collision time',
        'model':'noncrossing matching of exact published equal-time single-interval lengths',
        'algebra_checks':algebra, 'endpoint_cutoff':1e-4,
        'geometries':[x.tolist() for x in geometries],
        'radii':[.5,1.,2.], 'times':times.tolist(), 'party_counts':[3,4,5],
        'vectors':vectors, 'facet_evaluations':372*vectors,
        'facet_source':data['source'], 'facet_paper':data['paper'],
        'orbit_sizes':data['orbit_sizes'], 'minimum_per_orbit':orbit_min.tolist(),
        'minimum_facet_residual':all_min,'worst_case':worst_case,
        'minimum_discrete_kinematic_weight':weight_min,
        'kernel_reconstruction_max_error':reconstruction_error,
        'numerical_tolerance':1e-10,
        'example_curves':curves,
    }
    (args.output/'collision_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig, axs = plt.subplots(1,3,figsize=(12.2,3.7),layout='constrained')
    tt=np.array([r['t'] for r in curves])
    axs[0].plot(tt,[r['length_C_ren'] for r in curves],color='#16697a',lw=2.4)
    axs[0].set(title='Interval crossing an injection point',xlabel='BTZ boundary time $t$',ylabel='Renormalized geodesic length')
    axs[1].plot(tt,[r['minus_I3'] for r in curves],label='$-I_3(A:B:C)$',lw=2.2)
    axs[1].plot(tt,[r['cyclic5'] for r in curves],label='Five-party cyclic residual',lw=2.2)
    axs[1].axhline(0,color='.5',lw=.7)
    axs[1].set(title='Matching-model inequality residuals',xlabel='BTZ boundary time $t$',ylabel='Residual in length units')
    axs[1].legend(frameon=False,fontsize=8.5)
    theta=np.linspace(-PI,PI,360)
    axs[2].plot(np.cos(theta),np.sin(theta),color='.75',lw=1)
    colors=['#386cb0','#7f8c3a','#c0504d','#805eaa','#d58a24']
    for p in range(5):
        arc=np.linspace(example[2*p],example[2*p+1],90)
        axs[2].plot(np.cos(arc),np.sin(arc),lw=5,color=colors[p],solid_capstyle='round')
        mid=(example[2*p]+example[2*p+1])/2
        axs[2].text(1.18*math.cos(mid),1.18*math.sin(mid),'ABCDE'[p],ha='center',va='center',color=colors[p],weight='bold')
    axs[2].scatter([1,-1],[0,0],color='black',marker='*',s=75,zorder=5)
    axs[2].text(0,-1.47,'Stars: particle injection points\nAll 372 facets pass on the sampled grid',ha='center',fontsize=9)
    axs[2].set(title='Five boundary intervals',aspect='equal',xlim=(-1.4,1.4),ylim=(-1.6,1.4))
    axs[2].axis('off')
    fig.suptitle('Two-particle black-hole formation: exact interval kernel, candidate union matching',fontsize=12)
    fig.savefig(args.output/'collision_matching.png',dpi=190)
    fig.savefig(args.output/'collision_matching.pdf')
    print(json.dumps({k:v for k,v in report.items() if k not in ('times','geometries','example_curves')},indent=2))


if __name__=='__main__':
    main()
