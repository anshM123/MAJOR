/*
 * Reorganisation of the radial expansion of a d = 3 scalar block into the "lightcone" (rho, rhobar) monomial basis:
 *     Gt(r, eta) := sum_{n,j} B_{n,j} r^n P_j(eta) = (rho rhob)^{-J/2} sum_{i,q >= 0} C_{i,q} rho^i rhob^q ,
 *     rho = r e^{-chi}, rhob = r e^{chi}, eta = cosh chi,
 * using  P_j(cos th) = sum_{k=0}^{j} a_k a_{j-k} e^{i(j-2k)th},  a_k = (1/2)_k/k!.   With d = q - i, n = i + q - J:
 *     C_{i,q} = sum_{k'=0}^{min(i,q)} a_{k'} a_{|d|+k'} B_{n, |d|+2k'}     (terms with |j - J| > n absent).
 * Bt layout as in ads4/code/e2e/cblock.c:  B_{n,j} = Bt[off[j] + (n - |j-J|)/2]; it is first transposed to a level-major
 * array (B_{n,j} contiguous in j at fixed n, stride 2) for cache-friendly access.
 *     gcc -O3 -fopenmp -shared -fPIC lcconv.c -o liblcconv.so
 */
#include <stdlib.h>
#include <math.h>
/* Generalisation (derivative order l, lam = l + 1/2):  sum_{n,j} B_{n,j} r^n P_j^{(l)}(eta) = (2l-1)!! (rho rhob)^{-(J-l)/2} sum C^{(l)}_{i,q} rho^i rhob^q,
 * C^{(l)}_{i,q} = sum_{k'} a^lam_{k'} a^lam_{|d|+k'} B_{n, |d|+2k'+l},  n = i + q - J + l,  a^lam_k = (lam)_k/k!   (Gegenbauer:
 * P_j^{(l)} = (2l-1)!! C^{(l+1/2)}_{j-l},  C^{lam}_m(cos th) = sum_k a^lam_k a^lam_{m-k} e^{i(m-2k)th});  all terms non-negative.  */
void lc_convert_l(const double *Bt, const long *off, int J, int N, int I, int Q, int lsh, double *out);
void lc_convert(const double *Bt, const long *off, int J, int N, int I, int Q, double *out /* [(I+1)*(Q+1)] */)
{
    lc_convert_l(Bt, off, J, N, I, Q, 0, out);
}
void lc_convert_l(const double *Bt, const long *off, int J, int N, int I, int Q, int lsh, double *out)
{
    /* level-major copy: level n holds j = jlo(n), jlo+2, ..., J+n,  jlo = |J-n| if n <= J else (n-J)%2 ... use j0 = (J+n)%2 */
    long *loff = (long *)malloc((size_t)(N + 2)*sizeof(long));
    long tot = 0;
    for (int n = 0; n <= N; n++) { loff[n] = tot; tot += (J + n)/2 + 1; }      /* slots for j = (J+n)%2, ..., J+n step 2 */
    double *Lv = (double *)calloc((size_t)tot, sizeof(double));
    #pragma omp parallel for schedule(static)
    for (int n = 0; n <= N; n++) {
        int j0 = (J + n) % 2;
        for (int j = j0; j <= J + n; j += 2) {
            if (abs(j - J) > n) continue;
            Lv[loff[n] + (j - j0)/2] = Bt[off[j] + (n - abs(j - J))/2];
        }
    }
    int kmax = (I > Q ? I : Q) + J + N + 2;
    double *a = (double *)malloc((size_t)(kmax + 1)*sizeof(double));
    a[0] = 1.0;
    double lam = (double)lsh + 0.5;
    for (int k = 1; k <= kmax; k++) a[k] = a[k-1]*((double)(k - 1) + lam)/(double)k;
    #pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i <= I; i++) {
        for (int q = 0; q <= Q; q++) {
            int n = i + q - J + lsh;
            double s = 0.0;
            if (n >= 0 && n <= N) {
                int d = abs(q - i);
                int km = i < q ? i : q;
                int j0 = (J + n) % 2;                   /* j = d + 2k' + lsh has the parity of J + n */
                if ((d + lsh - j0) % 2 == 0) {
                    const double *L = Lv + loff[n] + (d + lsh - j0)/2;
                    /* |j - J| <= n:  j >= J - n  and  j <= J + n */
                    int kstart = 0;
                    if (d + lsh < J - n) kstart = (J - n - d - lsh + 1)/2;
                    int kend = km;
                    if (d + 2*kend + lsh > J + n) kend = (J + n - d - lsh)/2;
                    const double *ad = a + d;
                    for (int kp = kstart; kp <= kend; kp++) s += a[kp]*ad[kp]*L[kp];
                }
            }
            out[(long)i*(Q + 1) + q] = s;
        }
    }
    free(a); free(Lv); free(loff);
}
