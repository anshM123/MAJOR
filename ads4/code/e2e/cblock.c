/*
 * d = 3 scalar conformal blocks with general external dimensions: radial Casimir recursion and evaluation
 * of the radial series with (r d_r)- and eta-derivatives.  Generic in the real type:
 *     gcc -O3 -fopenmp -shared -fPIC cblock.c -o libcblock_d.so                 (double)
 *     gcc -O2 -fopenmp -shared -fPIC -DQUAD cblock.c -o libcblock_q.so -lquadmath   (__float128)
 * Reals cross the interface as (hi, lo) pairs of doubles (value = hi + lo).
 *
 *   g(r, eta) = (4 r)^Delta  sum_{n >= 0} sum_j  B_{n,j} r^n P_j(eta),   B_{0,J} = J!/(1/2)_J,
 * B_{n,j} != 0 only for |j - J| <= n, n + j - J even.  Storage: Bt[off[j] + i] = B_{|j-J|+2i, j}.
 * The recursion table (from casimir_terms.py) is a list of terms  r^k eta^i * coef * {1 | (1-eta^2) d_eta},
 * coef = c0 + c1 E + c2 E^2 + c3 lam + c4 C acting on r^E P_j (lam = j(j+1)).
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef QUAD
#include <quadmath.h>
typedef __float128 real;
#define EXP expq
#define LOG logq
#define COSH coshq
#define FABS fabsq
#else
typedef double real;
#define EXP exp
#define LOG log
#define COSH cosh
#define FABS fabs
#endif

static inline real from2(const double *p) { return (real)p[0] + (real)p[1]; }
static inline void to2(real x, double *p) { p[0] = (double)x; p[1] = (double)(x - (real)p[0]); }

long block_size(int J, int N) {
    long tot = 0;
    for (int j = 0; j <= J + N; j++) { int n0 = abs(j - J); if (n0 <= N) tot += (N - n0)/2 + 1; }
    return tot;
}

/* returns pointer to coefficient storage (to be freed with block_free) */
real *block_build(const double *Delta2, int J, int N, int nterms, const int *tk, const int *ti, const int *ttype,
                  const double *tc, long *off)
{
    real Delta = from2(Delta2);
    int jmax = J + N, nj = jmax + 1, pad = 8, L = nj + 2*pad;
    long tot = 0;
    for (int j = 0; j <= jmax; j++) { off[j] = tot; int n0 = abs(j - J); tot += (N - n0)/2 + 1; }
    real *Bt = (real *)calloc(tot, sizeof(real));
    real C = Delta*(Delta - 3) + (real)J*(J + 1);
    int depth = 9;
    real **lev = (real **)malloc(depth*sizeof(real *));
    for (int q = 0; q < depth; q++) lev[q] = (real *)calloc(L, sizeof(real));
    real *v = (real *)calloc(L, sizeof(real)), *w = (real *)calloc(L, sizeof(real)), *rhs = (real *)calloc(L, sizeof(real));
    real *eup = (real *)calloc(L, sizeof(real)), *edn = (real *)calloc(L, sizeof(real));
    real *dup = (real *)calloc(L, sizeof(real)), *ddn = (real *)calloc(L, sizeof(real));
    /* (eta v)_j = v_{j-1} j/(2j-1) + v_{j+1} (j+1)/(2j+3);  (D v)_j = v_{j+1} (j+1)(j+2)/(2j+3) - v_{j-1} (j-1) j/(2j-1) */
    for (int j = 0; j <= jmax + 4; j++) {
        real rj = (real)j;
        eup[j + pad] = (j >= 1) ? rj/(2*rj - 1) : 0;          /* multiplies v_{j-1} */
        edn[j + pad] = (rj + 1)/(2*rj + 3);                     /* multiplies v_{j+1} */
        dup[j + pad] = (j >= 1) ? -(rj - 1)*rj/(2*rj - 1) : 0;  /* multiplies v_{j-1} */
        ddn[j + pad] = (rj + 1)*(rj + 2)/(2*rj + 3);            /* multiplies v_{j+1} */
    }
    /* level 0 */
    real b0 = 1;
    for (int q = 0; q < J; q++) b0 *= (real)(q + 1)/((real)q + (real)0.5);
    lev[0][J + pad] = b0;
    Bt[off[J]] = b0;
    for (int n = 1; n <= N; n++) {
        int jlo = J - n; if (jlo < 0) jlo = 0;
        int jhi = J + n;
        memset(rhs, 0, L*sizeof(real));
        for (int tt = 0; tt < nterms; tt++) {
            int k = tk[tt];
            if (k == 0 || k > n) continue;
            real E = Delta + (real)(n - k);
            const double *c = tc + 5*tt;
            real *src = lev[(n - k) % depth];
            int slo = J - (n - k); if (slo < 0) slo = 0;
            int shi = J + (n - k);
            int lo = slo - 4; if (lo < 0) lo = 0;
            int hi = shi + 4; if (hi > jmax) hi = jmax;
            real base = (real)c[0] + (real)c[1]*E + (real)c[2]*E*E + (real)c[4]*C;
            for (int j = lo - 1 >= 0 ? lo - 1 : 0; j <= hi + 1 && j < nj + pad - 1; j++) v[j + pad] = 0;
            for (int j = slo; j <= shi; j++) v[j + pad] = (base + (real)c[3]*(real)j*(real)(j + 1))*src[j + pad];
            if (ttype[tt] == 1) {
                for (int j = lo; j <= hi; j++) w[j + pad] = ddn[j + pad]*v[j + 1 + pad] + dup[j + pad]*(j >= 1 ? v[j - 1 + pad] : 0);
                for (int j = lo; j <= hi; j++) v[j + pad] = w[j + pad];
                v[hi + 1 + pad] = 0; if (lo >= 1) v[lo - 1 + pad] = 0;
            }
            for (int q = 0; q < ti[tt]; q++) {
                for (int j = lo; j <= hi; j++) w[j + pad] = eup[j + pad]*(j >= 1 ? v[j - 1 + pad] : 0) + edn[j + pad]*v[j + 1 + pad];
                for (int j = lo; j <= hi; j++) v[j + pad] = w[j + pad];
                v[hi + 1 + pad] = 0; if (lo >= 1) v[lo - 1 + pad] = 0;
            }
            for (int j = lo; j <= hi; j++) rhs[j + pad] += v[j + pad];
        }
        real *dst = lev[n % depth];
        memset(dst, 0, L*sizeof(real));
        real E = Delta + (real)n;
        for (int j = jlo; j <= jhi; j++) {
            if (((n + j - J) & 1) != 0) continue;
            real K = E*(E - 3) + (real)j*(real)(j + 1) - C;
            real val = -rhs[j + pad]/(K/2);
            dst[j + pad] = val;
            Bt[off[j] + (n - abs(j - J))/2] = val;
        }
    }
    for (int q = 0; q < depth; q++) free(lev[q]);
    free(lev); free(v); free(w); free(rhs); free(eup); free(edn); free(dup); free(ddn);
    return Bt;
}

void block_free(real *Bt) { free(Bt); }

void block_get(const real *Bt, long idx, double *out2) { to2(Bt[idx], out2); }

/*
 * out[p][k][l] (as hi,lo pairs) = sum_{n,j} B_{n,j} sig^n (Delta+n)^k r^n P_j^{(l)}(cosh chi),  k <= kmax, l <= lmax,
 * sig = +1 or -1.  r, chi given as (hi, lo) pairs.  cut, pwcut: levels n with (r e^chi)^n e^{J chi} n^pwcut < exp(-cut) are dropped.
 */
void block_eval(const real *Bt, const long *off, const double *Delta2, int J, int N, int npts,
                const double *r2, const double *chi2, int sig, int kmax, int lmax, double cut, double pwcut, double *out)
{
    real Delta = from2(Delta2);
    int jmax = J + N;
    int nk = kmax + 1, nl = lmax + 1;
    #pragma omp parallel
    {
        real *Pt = (real *)malloc((size_t)(jmax + 2)*nl*sizeof(real));
        #pragma omp for schedule(dynamic, 4)
        for (int p = 0; p < npts; p++) {
            real r = from2(r2 + 2*p), chi = from2(chi2 + 2*p);
            real lr = LOG(r), emc = EXP(-chi), eta = COSH(chi), r2v = r*r;
            int nmaxp = N;
            /* terms decay like (r e^{chi})^n = rho_big^n */
            if (lr + chi < 0) {
                /* levels with (r e^chi)^n e^{J chi} n^pwcut < exp(-cut) are dropped; the observed truncation error of the conserved
                   tensors is e^{-X} X^p, X = n(-log rho_big), p ~ 5 (l = 1), so pwcut = 6..10 with cut = 30 is safe (t4b_convergence.py) */
                double pw = pwcut, lrb = (double)(lr + chi), base = cut + (double)J*(double)chi;
                double q = -base/lrb;
                for (int it = 0; it < 4; it++) q = -(base + pw*log(q + 10.0))/lrb;
                if (q < (double)N) nmaxp = (int)q + 1;
            }
            int jtop = J + nmaxp; if (jtop > jmax) jtop = jmax;
            /* scaled Legendre derivatives Pt[l][j] = P_j^{(l)}(eta) e^{-j chi} */
            for (int l = 0; l < nl; l++) { Pt[l*(jmax + 2) + 0] = (l == 0) ? 1 : 0; Pt[l*(jmax + 2) + 1] = (l == 0) ? eta*emc : ((l == 1) ? emc : 0); }
            for (int j = 1; j < jtop; j++) {
                Pt[j + 1] = ((2*(real)j + 1)*eta*emc*Pt[j] - (real)j*emc*emc*Pt[j - 1])/((real)j + 1);
                for (int l = 1; l < nl; l++)
                    Pt[l*(jmax + 2) + j + 1] = emc*emc*Pt[l*(jmax + 2) + j - 1] + (2*(real)j + 1)*emc*Pt[(l - 1)*(jmax + 2) + j];
            }
            real acc[8][8];
            for (int k = 0; k < nk; k++) for (int l = 0; l < nl; l++) acc[k][l] = 0;
            for (int j = 0; j <= jtop; j++) {
                int n0 = abs(j - J);
                if (n0 > nmaxp) continue;
                int imax = (nmaxp - n0)/2;
                const real *B = Bt + off[j];
                real w = EXP((real)n0*lr + (real)j*chi);
                if (sig < 0 && (n0 & 1)) w = -w;
                real s[8] = {0, 0, 0, 0, 0, 0, 0, 0};
                real x = Delta + (real)n0;
                for (int i = 0; i <= imax; i++) {
                    real t = B[i]*w;
                    real pw = t;
                    for (int k = 0; k < nk; k++) { s[k] += pw; pw *= x; }
                    w *= r2v; x += 2;
                }
                for (int k = 0; k < nk; k++) for (int l = 0; l < nl; l++) acc[k][l] += s[k]*Pt[l*(jmax + 2) + j];
            }
            for (int k = 0; k < nk; k++) for (int l = 0; l < nl; l++) to2(acc[k][l], out + 2*(((long)p*nk + k)*nl + l));
        }
        free(Pt);
    }
}
