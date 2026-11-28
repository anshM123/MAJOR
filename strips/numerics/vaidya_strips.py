"""Thin-null-shell Vaidya-AdS4 strips and all five-party RT cone facets.

Units: AdS radius L=1, final horizon zh=1; results are A/Ly=4GN*S/Ly.
The integrand is z^-2 sqrt(1-f*v'^2-2*v'*z'), not a 3D geodesic length.
Primary equations: Liu--Suh, arXiv:1311.1200, equations (3.3)--(3.23).
See NUMERICS.md for the domain of the branch exclusion and validation limits.
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar
from scipy.special import gamma
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parent
C = math.sqrt(math.pi) * gamma(0.75) / gamma(0.25)
EPSILON = 1e-4
NQUAD = 160


@lru_cache(None)
def gauss(n):
    x, w = leggauss(n)
    return (x + 1) / 2, w / 2


def moments(zt, r, n=NQUAD):
    """Return width/(2zt)-C, t/zt, and finite area correction*zt/2.

    u=z/zt=sqrt(sin(theta)), theta=asin(r^2)*s^2 desingularizes both ends.
    Subtraction uses (1/sqrt(Q)-1)=(1-Q)/(sqrt(Q)*(1+sqrt(Q))).
    """
    if r <= 0:
        return np.zeros(3)
    s, weights = gauss(n)
    theta_c = math.asin(min(1.0, r * r))
    # Resolve the boundary layer when the shell approaches the turning point.
    qc = math.sqrt(max(0.0, 1 - r**4))
    if 0 < qc < 0.05:
        scale = max(qc / (2 * theta_c), 1e-14)
        distances = [scale]
        while distances[-1] < 0.2:
            distances.append(min(0.2, distances[-1] * 8))
        bounds = sorted(set([0.0, 1.0] + [1 - x for x in distances]))
        ss, ww = [], []
        for lo, hi in zip(bounds[:-1], bounds[1:]):
            ss.append(lo + (hi - lo) * s)
            ww.append((hi - lo) * weights)
        s, weights = np.concatenate(ss), np.concatenate(ww)
    theta = theta_c * s * s
    u = np.sqrt(np.sin(theta))
    cos_theta = np.cos(theta)
    a = zt**3
    E = -0.5 * a * r * qc
    h = 1 - a * u**3
    energy_term = (E * u * u / cos_theta)**2
    Q = h + energy_term
    root = np.sqrt(Q)
    ratio = (a * u**3 - energy_term) / (root * (1 + root))
    measure = weights * theta_c * s
    return np.array([
        np.sum(measure * u * ratio),
        np.sum(measure * (cos_theta / u + E * u / root) / h),
        np.sum(measure * ratio / u**3),
    ])


def finite_area(zt, r, correction, epsilon=EPSILON):
    """Exact radial truncation to Gaussian accuracy, removing only 2/epsilon.

    Both the vacuum small-u integral and the shell correction below the cutoff
    are subtracted. Thus this is not an arbitrary vacuum-subtracted cone vector.
    """
    s, weights = gauss(24)
    delta = epsilon / zt
    u = delta * s
    B = 1 - u**4
    sqB = np.sqrt(B)
    vacuum_lower = delta * np.dot(weights, u*u / (sqB * (1 + sqB)))
    upper = min(delta, r)
    if upper:
        u = upper * s
        a = zt**3
        E2 = 0.25 * a*a * r*r * max(0.0, 1 - r**4)
        B = 1 - u**4
        F = (1 - a*u**3)*B + E2*u**4
        sqB, sqF = np.sqrt(B), np.sqrt(F)
        delta_integrand = (a*u*B - E2*u*u) / (sqF*sqB*(sqF + sqB))
        shell_lower = upper * np.dot(weights, delta_integrand)
    else:
        shell_lower = 0.0
    return 2 / zt * (-C + correction - vacuum_lower - shell_lower)


class Strip:
    def __init__(self, width, n=NQUAD):
        self.width, self.n = float(width), n
        self.zvac = self.width / (2*C)
        self.calls = 0
        self.curve = lru_cache(maxsize=None)(self._curve)
        self.zstatic, self.tsat, self.astatic = self.curve(1.0)
        self.avac = finite_area(self.zvac, 0, 0)
        # Cover the whole exterior family, including the near-saturation cusp.
        rr = np.unique(np.round(np.r_[np.linspace(0, 0.98, 99),
                             1 - np.geomspace(1e-14, .02, 100), 1.0],15))
        self.rscan = np.sort(rr)
        self.scan = np.array([self.curve(float(r)) for r in self.rscan])
        self.tscan = self.scan[:, 1]
        peak=int(np.argmax(self.tscan))
        if 0<peak<len(self.rscan)-1 and self.tscan[peak]-self.tsat>1e-10:
            lo,hi=self.rscan[peak-1],self.rscan[peak+1]
            if hi<1:
                result=minimize_scalar(lambda p:-self.curve(1-10**p)[1],
                    bounds=(math.log10(1-hi),math.log10(1-lo)),method="bounded",
                    options={"xatol":2e-10})
                rpeak=1-10**result.x
                self.rscan=np.sort(np.r_[self.rscan,rpeak])
                self.scan=np.array([self.curve(float(r)) for r in self.rscan])
                self.tscan=self.scan[:,1]
        self.monotone = bool(np.all(np.diff(self.tscan) > 0))
        self.tmax = float(np.max(self.tscan))
        self.max_candidates = 1

    def _curve(self, r):
        self.calls += 1
        if r <= 0:
            return self.zvac, 0.0, finite_area(self.zvac, 0, 0)
        def equation(zt):
            dw = moments(zt, r, self.n)[0]
            return 2*zt*(C + dw) - self.width
        zt = brentq(equation, self.zvac*.25, self.zvac,
                    xtol=2e-14, rtol=2e-14)
        dw, tau, correction = moments(zt, r, self.n)
        return zt, zt*tau, finite_area(zt, r, correction)

    def at_time(self, t):
        if t <= 0:
            return self.avac, self.zvac, 0.0, "vacuum"
        candidates=[]
        if t >= self.tsat:
            candidates.append((self.astatic,self.zstatic,1.0,"static"))
        if t <= self.tmax:
            values=self.tscan-t
            crossings=np.flatnonzero(values[:-1]*values[1:]<0)
            roots=[]
            for i in crossings:
                lo,hi=float(self.rscan[i]),float(self.rscan[i+1])
                roots.append(brentq(lambda x:self.curve(x)[1]-t,lo,hi,
                                    xtol=2e-15,rtol=2e-15))
            roots.extend(float(self.rscan[i]) for i in np.flatnonzero(values==0))
            for r in roots:
                zt,_,area=self.curve(r)
                candidates.append((area,zt,r,"shell"))
        if not candidates:
            raise RuntimeError(f"No extremal branch found for w={self.width},t={t}")
        self.max_candidates=max(self.max_candidates,len(candidates))
        return min(candidates,key=lambda candidate:candidate[0])


def make_intervals(widths, gap=0.0):
    x, out = 0.0, []
    for width in widths:
        out.append((round(x, 10), round(x+width, 10)))
        x += width + gap
    return out


def merged_intervals(intervals, mask):
    selected = [intervals[i] for i in range(len(intervals)) if mask & (1 << i)]
    result = []
    for lo, hi in selected:
        if result and abs(result[-1][1]-lo) < 1e-11:
            result[-1] = (result[-1][0], hi)
        else:
            result.append((lo, hi))
    return result


def subset_mask(label):
    return sum(1 << (ord(c)-ord("A")) for c in label)


def union_area(intervals, mask, lookup):
    """Catalan interval DP. Each matching has m pairs, hence UV term 2m/eps."""
    merged = merged_intervals(intervals, mask)
    endpoints = tuple(x for interval in merged for x in interval)
    @lru_cache(None)
    def dp(i, j):
        if i > j:
            return 0.0
        return min(lookup[round(endpoints[k]-endpoints[i], 10)]
                   + dp(i+1, k-1) + dp(k+1, j)
                   for k in range(i+1, j+1, 2))
    return dp(0, len(endpoints)-1), len(merged)


def width_catalog(configs):
    widths = set()
    for intervals in configs.values():
        for mask in range(1, 1 << len(intervals)):
            ee = [x for pair in merged_intervals(intervals, mask) for x in pair]
            for i in range(len(ee)):
                for j in range(i+1, len(ee), 2):
                    widths.add(round(ee[j]-ee[i], 10))
    return sorted(widths)


def check_multiple_starts(strips):
    """Independent bounded 2D root search at representative times/widths."""
    checks = []
    selected = [strips[0], strips[len(strips)//2], strips[-1]]
    max_zt_difference, max_r_difference, failed = 0.0, 0.0, 0
    for strip in selected:
        for fraction in [.13, .47, .83, .97]:
            target = fraction*strip.tsat
            _, refzt, refr, _ = strip.at_time(target)
            roots = []
            for r0 in [.08, .35, .7, .95]:
                def residual(x):
                    zt, r = x
                    dw, tau, _ = moments(zt, r, NQUAD)
                    return [(2*zt*(C+dw)-strip.width)/strip.width,
                            (zt*tau-target)/strip.tsat]
                sol = least_squares(residual, [strip.zvac*.97, r0],
                                    bounds=([strip.zvac*.25, 1e-9],
                                            [strip.zvac*1.05, 1-1e-12]),
                                    xtol=5e-14, ftol=5e-14, gtol=5e-14,
                                    max_nfev=120)
                error = float(np.max(np.abs(residual(sol.x))))
                if error > 3e-8:
                    failed += 1
                roots.append({"initial_r":r0,"zt":float(sol.x[0]),
                              "r":float(sol.x[1]),"scaled_residual":error})
                max_zt_difference = max(max_zt_difference, abs(sol.x[0]-refzt))
                max_r_difference = max(max_r_difference, abs(sol.x[1]-refr))
            checks.append({"width":strip.width,"time":target,"roots":roots})
    return {"checks":checks,"failed_roots":failed,
            "max_zt_difference_from_continuation":max_zt_difference,
            "max_r_difference_from_continuation":max_r_difference}


def main():
    started = time.time()
    facets_data = json.loads((ROOT/"data"/"facets_n5_all372.json").read_text())
    facets = np.array(facets_data["facets"], dtype=np.int64)
    labels = facets_data["subset_order"]
    contiguous = make_intervals([.10,.14,.18,.16,.12])
    disjoint = make_intervals([.12]*5, .0741)
    configs = {"contiguous":contiguous,"disjoint":disjoint}
    widths = width_catalog(configs)
    print(f"Building {len(widths)} distinct interval curves", flush=True)
    strips = [Strip(w) for w in widths]
    by_width = dict(zip(widths, strips))
    tsat_max = max(s.tsat for s in strips)
    times = np.linspace(-.025, tsat_max*1.06, 401)
    # Resolve each possible first-order saturation window as well as the dense grid.
    transition_times=np.concatenate([np.linspace(s.tsat,s.tmax,7)
                                    for s in strips if s.tmax-s.tsat>1e-9])
    times = np.unique(np.r_[times,0.0,[s.tsat for s in strips],transition_times])
    finite = np.empty((len(times),len(widths)))
    tips = np.empty_like(finite)
    ratios = np.empty_like(finite)
    regimes = np.empty_like(finite, dtype=np.int8)
    residual_w, residual_t = 0.0, 0.0
    for j, strip in enumerate(strips):
        for i, t in enumerate(times):
            a, zt, r, regime = strip.at_time(float(t))
            finite[i,j], tips[i,j], ratios[i,j] = a, zt, r
            regimes[i,j] = {"vacuum":0,"shell":1,"static":2}[regime]
            if regime == "shell":
                dw, tau, _ = moments(zt, r, NQUAD)
                residual_w = max(residual_w, abs(2*zt*(C+dw)-strip.width))
                residual_t = max(residual_t, abs(zt*tau-t))
        print(f"  width={strip.width:.4f} saturation={strip.tsat:.6f}", flush=True)
    summaries, arrays = {}, {}
    for config_name, full_intervals in configs.items():
        for n in [3,4,5]:
            intervals = full_intervals[:n]
            ent = np.zeros((len(times),31))
            components = np.zeros(31,dtype=np.int64)
            mask_cache = {}
            for k, label in enumerate(labels):
                mask = subset_mask(label) & ((1<<n)-1)
                if mask not in mask_cache:
                    values=[]
                    for i in range(len(times)):
                        lookup=dict(zip(widths,finite[i]))
                        a, m = union_area(intervals,mask,lookup)
                        values.append(a)
                    mask_cache[mask] = (np.array(values),m)
                ent[:,k],components[k] = mask_cache[mask]
            uv_coeff = facets@components
            assert np.min(uv_coeff)>=0
            finite_slacks = ent@facets.T
            slacks = finite_slacks + (2/EPSILON)*uv_coeff[np.newaxis,:]
            safe = uv_coeff==0
            min_all = float(np.min(slacks))
            min_finite = float(np.min(finite_slacks[:,safe]))
            violations = int(np.count_nonzero(slacks < -2e-8))
            key=f"{config_name}_{n}"
            orbit_minima=[]
            for orbit in facets_data["facet_orbits"]:
                finite_orbit=[j for j in orbit if safe[j]]
                orbit_minima.append(np.min(finite_slacks[:,finite_orbit],axis=1)
                                    if finite_orbit else np.full(len(times),np.nan))
            arrays[key+"_finite_entropies"]=ent
            arrays[key+"_components"]=components
            arrays[key+"_finite_facet_slacks"]=finite_slacks
            arrays[key+"_uv_coefficients"]=uv_coeff
            arrays[key+"_orbit_minima"]=np.array(orbit_minima).T
            arrays[key+"_facet_minima"]=np.min(slacks,axis=1)
            summaries[key]={"parties":n,"intervals":intervals,
                            "facet_rows_tested":372,
                            "nontrivial_projected_rows":int(np.count_nonzero(np.max(np.abs(slacks),axis=0)>1e-10)),
                            "UV_finite_rows":int(np.count_nonzero(safe)),
                            "minimum_regulated_slack":min_all,
                            "minimum_UV_finite_slack":min_finite,
                            "violations_below_minus_2e-8":violations}
    # Independent quadrature refinement at a spread of shell points.
    max_refinement_w=max_refinement_t=max_refinement_a=0.0
    refinement_count=0
    for j in np.unique(np.linspace(0,len(widths)-1,8,dtype=int)):
        strip=strips[j]
        for fraction in [.01,.15,.43,.79,.97,.999]:
            a, zt, r, _=strip.at_time(float(fraction*strip.tsat))
            dw,tau,correction=moments(zt,r,2*NQUAD)
            anew=finite_area(zt,r,correction)
            max_refinement_w=max(max_refinement_w,abs(2*zt*(C+dw)-strip.width))
            max_refinement_t=max(max_refinement_t,abs(zt*tau-fraction*strip.tsat))
            max_refinement_a=max(max_refinement_a,abs(anew-a))
            refinement_count+=1
    multistart=check_multiple_starts(strips)
    # Static saturation must agree with the v_tip=0 exterior static solution.
    from scipy.integrate import quad
    static_time_error=max(abs(s.tsat-quad(lambda z:1/(1-z**3),0,s.zstatic,
                                        epsabs=1e-12,epsrel=1e-12)[0]) for s in strips)
    validation={
        "model":"thin-null-shell planar Vaidya-AdS4",
        "metric":"ds2=z^-2[-(1-Theta(v)z^3)dv^2-2dv dz+dx^2+dy^2]",
        "area_normalization":"A/Ly=4GN*S/Ly, L=zh=1",
        "source":"https://arxiv.org/abs/1311.1200, section III.A",
        "facet_source":"https://arxiv.org/abs/1903.09148; author HEC database",
        "epsilon":EPSILON,"gauss_order":NQUAD,
        "time_samples":len(times),"time_min":float(times[0]),"time_max":float(times[-1]),
        "interval_width_count":len(widths),"max_width":max(widths),
        "max_vacuum_tip_bound":float(max(s.zvac for s in strips)),
        "minimum_exterior_h_lower_bound":float(1-max(s.zvac for s in strips)**3),
        "branch_selection":"All sampled t(r) intersections plus the admissible static branch; smallest regulated area selected.",
        "all_sampled_parameterized_branches_monotone":all(s.monotone for s in strips),
        "maximum_extremal_candidates_compared":max(s.max_candidates for s in strips),
        "saturation_windows":[{"width":s.width,"static_admissibility_time":s.tsat,
                               "maximum_shell_time":s.tmax,
                               "overshoot":s.tmax-s.tsat,
                               "maximum_candidates":s.max_candidates} for s in strips],
        "minimum_parameterized_time_increment":float(min(np.min(np.diff(s.tscan)) for s in strips)),
        "max_width_equation_residual":float(residual_w),
        "max_time_equation_residual":float(residual_t),
        "quadrature_refinement":{"points":refinement_count,"orders":[NQUAD,2*NQUAD],
            "max_width_difference":float(max_refinement_w),
            "max_time_difference":float(max_refinement_t),
            "max_finite_area_difference":float(max_refinement_a)},
        "static_time_independent_quad_error":float(static_time_error),
        "multiple_start_root_audit":multistart,
        "configuration_results":summaries,
        "interpretation":"Sampled floating-point consistency check; not continuous-time or interval-arithmetic certification.",
        "elapsed_seconds":time.time()-started,
    }
    connection=lambda t:(2*by_width[.12].at_time(t)[0]
                         -by_width[.0741].at_time(t)[0]
                         -by_width[.3141].at_time(t)[0])
    transition=brentq(connection,.10,.30,xtol=2e-13)
    validation["two_strip_mutual_information_transition"]={
        "width":.12,"gap":.0741,"time":transition,
        "signed_connection_margin_at_root":connection(transition),
        "vacuum_mutual_information":connection(0),
        "late_signed_connection_margin":connection(.5)}
    np.savez_compressed(ROOT/"vaidya_data.npz",times=times,widths=widths,
                        strip_finite_areas=finite,strip_tips=tips,shell_ratios=ratios,
                        regimes=regimes,facet_coefficients=facets,**arrays)
    (ROOT/"validation.json").write_text(json.dumps(validation,indent=2)+"\n")
    (ROOT/"configurations.json").write_text(json.dumps({"configs":configs,
            "subset_order":labels,"facet_orbit_sizes":facets_data["orbit_sizes"]},indent=2)+"\n")
    plot(times,widths,finite,arrays,facets_data,validation)
    print(json.dumps({"elapsed_seconds":validation["elapsed_seconds"],
                      "refinement":validation["quadrature_refinement"],
                      "results":summaries},indent=2),flush=True)


def plot(times,widths,finite,arrays,facet_data,validation):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":10,"axes.titlesize":11,"axes.labelsize":10,
                         "legend.fontsize":8,"figure.dpi":160,"savefig.dpi":220,
                         "font.family":"DejaVu Sans","axes.spines.top":False,
                         "axes.spines.right":False})
    fig, axes=plt.subplots(2,2,figsize=(10.0,7.2),constrained_layout=True)
    ax=axes[0,0]
    choices=[.12,.24,.42,.70,.8964]
    for w in choices:
        j=int(np.argmin(np.abs(np.array(widths)-w)))
        actual=widths[j]
        delta=finite[:,j]-finite[0,j]
        ax.plot(times,delta,label=f"$w={actual:.4g}$",lw=1.8)
    ax.set(title="(a) Strip entropy growth",xlabel="$t/z_h$",ylabel=r"$4G_N\,\Delta S/L_y$")
    ax.legend(frameon=False,ncol=2)
    ax.axvline(0,color=".6",ls=":",lw=.8)
    ax=axes[0,1]
    labels=facet_data["subset_order"]
    e=arrays["disjoint_5_finite_entropies"]
    mi=e[:,labels.index("A")]+e[:,labels.index("B")]-e[:,labels.index("AB")]
    ax.plot(times,mi,color="#c05d24",lw=2)
    ax.axhline(0,color=".5",lw=.8)
    ax.set(title="(b) Two-strip mutual information",xlabel="$t/z_h$",
           ylabel=r"$4G_N\,I(A:B)/L_y$")
    ax.text(.98,.95,"Five equal strips: width 0.12\ngap 0.0741",transform=ax.transAxes,
            ha="right",va="top",fontsize=9)
    ax=axes[1,0]
    # Display each orbit's time-dependent spread, since many actual minima vanish.
    ss=arrays["contiguous_5_finite_facet_slacks"]
    uv=arrays["contiguous_5_uv_coefficients"]
    for k,orbit in enumerate(facet_data["facet_orbits"]):
        jj=[j for j in orbit if uv[j]==0]
        if not jj:continue
        maximum=np.max(ss[:,jj],axis=1)
        minimum=np.min(ss[:,jj],axis=1)
        ax.plot(times,maximum,lw=1.4,label=f"orbit {k+1}")
        ax.plot(times,np.maximum(minimum,0),lw=.7,alpha=.45)
    ax.set(title="(c) UV-finite facet margins: orbit envelopes",
           xlabel="$t/z_h$",ylabel="Upper / lower margin in each orbit")
    ax.set_yscale("symlog",linthresh=1e-3)
    ax.legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,1.0))
    ax=axes[1,1]
    names=["contiguous_3","contiguous_4","contiguous_5",
           "disjoint_3","disjoint_4","disjoint_5"]
    for i,name in enumerate(names):
        margin=arrays[name+"_facet_minima"]
        ax.plot(times,margin/1e-12,label=name.replace("_"," / "),lw=1.1)
    ax.axhline(0,color=".35",lw=.8)
    ax.set(title="(d) Minimum of all 372 tested facet rows",
           xlabel="$t/z_h$",ylabel="Minimum margin $/10^{-12}$")
    ax.legend(frameon=False,ncol=2,fontsize=7)
    for ax in axes.flat:
        ax.grid(alpha=.18)
        ax.set_xlim(times[0],times[-1])
    fig.suptitle("Thin-shell thermalization in AdS$_4$: strip entropy cone check",fontsize=14)
    fig.savefig(ROOT/"vaidya_diagnostics.png")
    fig.savefig(ROOT/"vaidya_diagnostics.pdf")
    plt.close(fig)
    # A compact physics figure; the complete numerical audit remains separate.
    fig,axes=plt.subplots(1,2,figsize=(10.2,3.75),constrained_layout=True)
    ax=axes[0]
    for w in choices:
        j=int(np.argmin(np.abs(np.array(widths)-w)))
        ax.plot(times,finite[:,j]-finite[0,j],lw=1.9,
                label=f"$w={widths[j]:.4g}$")
    ax.set(title="(a) Strip entanglement after the quench",xlabel="$t/z_h$",
           ylabel=r"$4G_N\,[S(t)-S_{\rm vac}]/L_y$")
    ax.legend(frameon=False,ncol=2,loc="upper left")
    ax.axvline(0,color=".6",ls=":",lw=.8)
    ax=axes[1]
    ax.plot(times,mi,color="#c05d24",lw=2)
    ax.axhline(0,color=".5",lw=.8)
    ax.set(title="(b) Mutual information crosses a phase boundary",
           xlabel="$t/z_h$",ylabel=r"$4G_N\,I(A:B)/L_y$")
    ax.text(.98,.95,"Five strips: width 0.12, gap 0.0741",transform=ax.transAxes,
            ha="right",va="top",fontsize=9)
    for ax in axes:
        ax.grid(alpha=.18)
        ax.set_xlim(times[0],times[-1])
    fig.suptitle("AdS$_4$ thin-shell thermalization",fontsize=13)
    fig.savefig(ROOT/"vaidya_entropy_cone.png")
    fig.savefig(ROOT/"vaidya_entropy_cone.pdf")
    plt.close(fig)


if __name__=="__main__":
    main()
