/*
 * Reorganisation of the radial expansion of a d = 3 scalar block into the "lightcone" (rho, rhobar) monomial basis:
 *     Gt(r, eta) := sum_{n,j} B_{n,j} r^n P_j(eta) = (rho rhob)^{-J/2} sum_{i,q >= 0} C_{i,q} rho^i rhob^q ,
 *     rho = r e^{-chi}, rhob = r e^{chi}, eta = cosh chi,
 * using  P_j(cos th) = sum_{k=0}^{j} a_k a_{j-k} e^{i(j-2k)th},  a_k = (1/2)_k/k!.   With d = q - i, n = i + q - J:
 *     C_{i,q} = sum_{k'=0}^{min(i,q)} a_{k'} a_{|d|+k'} B_{n, |d|+2k'}     (terms with |j - J| > n absent).
 * Bt layout as in ads4/code/e2e/cblock.c:  B_{n,j} = Bt[off[j] + (n - |j-J|)/2].
 *     gcc -O3 -fopenmp -shared -fPIC lcconv.c -o liblcconv.so
 */
#include <stdlib.h>
#include <math.h>
void lc_convert(const double *Bt, const long *off, int J, int N, int I, int Q, double *out /* [(I+1)*(Q+1)] */)
{
    int kmax = (I > Q ? I : Q) + J + N + 2;
    double *a = (double *)malloc((size_t)(kmax + 1)*sizeof(double));
    a[0] = 1.0;
    for (int k = 1; k <= kmax; k++) a[k] = a[k-1]*((double)k - 0.5)/(double)k;
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i <= I; i++) {
        for (int q = 0; q <= Q; q++) {
            int n = i + q - J;
            double s = 0.0;
            if (n >= 0 && n <= N) {
                int d = abs(q - i);
                int km = i < q ? i : q;
                for (int kp = 0; kp <= km; kp++) {
                    int j = d + 2*kp;
                    if (abs(j - J) > n) continue;
                    s += a[kp]*a[d + kp]*Bt[off[j] + (n - abs(j - J))/2];
                }
            }
            out[(long)i*(Q + 1) + q] = s;
        }
    }
    free(a);
}
