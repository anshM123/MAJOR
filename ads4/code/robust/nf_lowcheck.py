"""Exact low side of a no-forward functional recomputed from the amplitude (2.7a) and matter (B.1) by residues: lowcheck.py with the
coefficients taken from functionals_noforward.json (via nf_cert.load).  Usage: python3 nf_lowcheck.py NAME"""
import sys
import rcert, nf_cert
rcert.load_functional = nf_cert.load
import lowcheck
lowcheck.main(sys.argv[1])
