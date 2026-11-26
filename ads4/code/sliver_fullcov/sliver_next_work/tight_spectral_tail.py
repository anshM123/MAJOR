from fractions import Fraction as F
import csv,sys,json
from base_bounds import logl,logu,sqrtlo,sqrtup
from resid_budget import pow2exp,gaussian_mass
# Use retained absolute residual panel masses (psi2 is supported on [0,Lambda]).
# Add a conservative low-frequency bound and the two cutoff corrections separately.
# For this proof the weight g = ideal + residual, residual below 1040.1 from existing ledger.
# The unweighted cutoff correction is only known globally; its worst frequency may be 1040.1.
B=F(11,7);R=F(199,200);H=F(1,998**2)
rows=list(csv.DictReader(open('band_masses.csv')))
weighted=F(120000001)*pow2exp(B*999)+F(60000030)*pow2exp(B*160)
weighted+=F(2,998)*sum(F(x['mass'])*pow2exp(B*F(x['hi'])) for x in rows)
weighted+=F(2,998)*F('1.34e-14')*pow2exp(B*F('1040.1'))
weighted+=F(2,998)*gaussian_mass(F('1040'),B)
M=2*(1-R)**-4*weighted
# bound by exact power of 2
E=M.numerator.bit_length()-M.denominator.bit_length()+1
while M>2**E:E+=1
assert M<=2**E
print('exponent',E,'old 2417')
def tail(c,K):
 t=F(c)/R;assert 0<t<1
 return 2**E*t**(K+1)*((1+20*H)/(1-t)+2*H*((K+5)**2/(1-t)+2*(K+5)*t/(1-t)**2+t*(1+t)/(1-t)**3))
