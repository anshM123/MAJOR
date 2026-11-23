"""Target 5: conserved <T T O_{Delta,l}> structure counts in d = 3 for generic Delta, by helicity counting
(conserved spin-2 operators carry helicities h = +-2; a generic spin-l operator O admits the pair (h1, h2) iff |h1 - h2| <= l;
parity maps (h1, h2) -> (-h1, -h2); Bose symmetry 1 <-> 2 maps (h1, h2) -> (h2, h1) with sign (-1)^l).
Each orbit of the pair set under parity gives one structure per parity eigenvalue that survives the Bose projection.
Parity label convention: in the (2,2) orbit the two eigen-combinations are one even and one odd; in the (2,-2) orbit the surviving
combination has Bose eigenvalue s(-1)^l = +1 and parity eigenvalue s, s = (-1)^l, i.e. parity (-1)^l (fixed by l = 4: DKKPS 2+ + 1-).
Output: cache/t5_counts.txt"""
out = open('cache/t5_counts.txt', 'w')
def log(s): print(s); out.write(s + '\n')
def count(l):
    even = odd = 0
    # orbit {(2,2), (-2,-2)}: Bose maps each pair to itself with sign (-1)^l
    if l % 2 == 0:          # both parity combinations survive
        even += 1; odd += 1
    # orbit {(2,-2), (-2,2)}: needs |h1 - h2| = 4 <= l
    if l >= 4:
        s = (-1)**l          # surviving combination (2,-2) + s(-2,2)
        if s == 1: even += 1
        else: odd += 1
    return even, odd
log('l : even + odd   (DKKPS 1708.05718 table: O_0 1+ + 1-, O_2 1+ + 1-, O_{2n>=4} 2+ + 1-, O_{2n+1>=5} 1-; P2_SPINNING §1.3: 1+1, 0, 1+1, 0, 2+1, 0+1, 2+1, 0+1, 2+1, 0+1)')
for l in range(0, 12):
    e, o = count(l)
    log(f'{l:2d}: {e} + {o}')
