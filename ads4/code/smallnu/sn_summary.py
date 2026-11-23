"""
sn_summary (SMALLNU_WINDOW.md): one line per state from cache/sn_odd_J*.json and cache/sn_even_J*.json -> cache/sn_summary.txt
"""
import os, sys, json, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import sn_guard  # noqa: F401
C = os.path.join(HERE, 'cache')
L = []
for f in sorted(glob.glob(os.path.join(C, 'sn_odd_J*.json')), key=lambda s: int(s.split('_J')[1][:-5])):
    d = json.load(open(f))
    L.append(f"J={d['J']} (odd, row {d['D']:g}): nu_min {d['numin']:.1f}; window max(|G|+3s) {d['G_meas']:.0f}, max|G| {d['G_cent']:.0f}")
    for k, o in d['states'].items():
        L.append(f"   m/L={o['mfl']:.4f} lam_c {o['lamc']:+.3f} (signed {o['lam_signed']:+.3f}) | base budget {o['E_base']:.3f} (trunc {o['E_trunc']:.3f}) | "
                 f"window conv {o['E_conv']:.3f} env1 {o['E_env1']:.3f} env3 {o['E_env3']:.3f} env10 {o['E_env10']:.3f} signed {o['E_lo'] + o['E_signed']:.3f} "
                 f"(<nu_min {o['E_lo']:.3f}) | L conv {o['L_conv']:.2f} env3 {o['L_env3']:.2f} env10 {o['L_env10']:.2f} signed {o['L_signed']:.2f} | "
                 f"cal-trunc L env3 {o['Lcal_env3']:.2f} signed {o['Lcal_signed']:.2f}")
for f in sorted(glob.glob(os.path.join(C, 'sn_even_J*.json')), key=lambda s: int(s.split('_J')[1][:-5])):
    d = json.load(open(f))
    L.append(f"J={d['J']} (even, rows {d['rows']}): nu_min {d['numin']:.1f}; window diag max(|G|+3s) {d['G_meas']:.0f}, max|G| {d['G_cent']:.0f}, "
             f"off {d['G_meas_off']:.1f}")
    for k, o in d['states'].items():
        s = (f"   m/L={o['mfl']:.4f} window conv {o['E_conv']:.3f} env1 {o['E_env1']:.3f} env3 {o['E_env3']:.3f} env10 {o['E_env10']:.3f} "
             f"measured {o['E_meas']:.3f} (<nu_min {o['E_lo']:.3f}; |dF|2 {o['dFwin_spec']:.4f}, err {o['E_signed']:.4f}; first-order {o['dlam_win_hyb']:+.4f})")
        if 'lamc' in o:
            s += (f" | lam_c {o['lamc']:+.3f} other {o['rest']:.3f} (lam-level {o['rest_lam']:.3f}) | L old {o['L_old']:.2f} -> env3 "
                  f"{o['L_env3']:.2f} env10 {o['L_env10']:.2f} meas {o['L_meas']:.2f}; lam-level old {o['L_old_lam']:.2f} -> "
                  f"{o['Llam_env3']:.2f} / {o['Llam_env10']:.2f} / {o['Llam_meas']:.2f}")
        L.append(s)
print("\n".join(L))
open(os.path.join(C, 'sn_summary.txt'), 'w').write("\n".join(L) + "\n")
