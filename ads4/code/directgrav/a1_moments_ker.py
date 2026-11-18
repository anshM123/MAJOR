"""corrected tensor-kernel circle weights (copy of the table in a1_moments.py)"""
KER = {
    'E': {(0, 0): 1/72, (2, -2): 0.25, (-2, 2): 0.25},                                            # D1^2 (Omega^U even kernel)
    'FF': {(0, 0): 1/72, (2, -2): 0.25, (-2, 2): 0.25, (1, -1): -1/32, (-1, 1): -1/32},           # S_ff = D1^2 + D1 D2 (flip-free)
    'KS': {(2, -2): 0.25, (-2, 2): 0.25, (1, -1): -1/64, (-1, 1): -1/64},                         # K* = (D1 + D2)(D1 - D2/2)
    'O': {(2, -2): 1/(4j), (-2, 2): -1/(4j)},                                                     # A = D5 D1 (odd kernel)
}


