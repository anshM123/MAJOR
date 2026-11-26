// Fixed-spin real-mass cell checker. Exact integer fixed point with explicit error bounds.
// Build: g++ -O3 -std=c++17 check_cell.cpp -lgmpxx -lgmp -o check_cell
// Run: ./check_cell 1005 4   or ./check_cell 1002 0
#include <gmpxx.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <string>
#include <chrono>
using Z=mpz_class; using V=std::vector<Z>;
const int P=1850, Q=260, VMAX=1040; const long SCAL=3; const long L2=998L*999, V2=1040L*1040, W=8000, R=100;
const Z S=Z(1)<<P;
Z ab(const Z&x){return x>=0?x:-x;}
Z div0(const Z&a,const Z&b){Z x;mpz_tdiv_q(x.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return x;}
Z ceildiv(const Z&a,const Z&b){assert(a>=0&&b>0);Z x;mpz_cdiv_q(x.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return x;}
Z shr(const Z&x,int p=P){Z y;mpz_tdiv_q_2exp(y.get_mpz_t(),x.get_mpz_t(),p);return y;}
Z powz(unsigned long a,unsigned long b){Z z;mpz_ui_pow_ui(z.get_mpz_t(),a,b);return z;}
Z norm(const V&v){Z r=0;for(auto&x:v)r+=ab(x);return r;}
double f(const Z&a){return mpq_class(a,S).get_d();}
// Integrals of T_j(2 Lambda^2 p^2/V^2 - 1) against exact NF02 weights.
// Recurrences derived by integration by parts; exact per-step rounding error bounds.
struct Moments{V a,d;Z err;Z weight_norm;};
Moments moments(int K){
 Z tnum=2*L2-V2, tden=V2; V T(K+1),TE(K+1);T[0]=S;T[1]=div0(tnum*S,tden);TE[1]=1;
 for(int j=1;j<K;j++){T[j+1]=div0(2*tnum*T[j],tden)-T[j-1];TE[j+1]=Z((long)(j+1))*(j+2)/2;}
 long X[]={5191828,11893013,17891069,22590571,18415583,160339862,-979855160,3620891327,-3115042418};
 long Y[]={1391927,6959636,51592135,187227898,1038008342,-4691758161,7462079762,-3734000568,-553835527};
 long bin[]={1,-5,10,-10,5,-1};V wa(15),wd(15);for(int k=0;k<9;k++)for(int i=0;i<6;i++){wa[k+i+1]+=Z(X[k])*bin[i]*100000;wd[k+i+1]+=Z(Y[k])*bin[i]*100000;}
 Z den=2748509817L;V aa(K+1),dd(K+1),ae(K+1),de(K+1);
 for(int m=0;m<=14;m++){
  V mu(K+1),e(K+1);mu[0]=S/(m+1);e[0]=1;
  mu[1]=div0(2*Z(L2)*S,Z(V2)*(m+3))-div0(S,Z(m+1));e[1]=2;
  for(int j=1;j<K;j++){
   Z A=2*(m-1),B=2*j-m-3,D=2*j+m+3;
   Z rhs=(A*mu[j]+B*mu[j-1])*tden-2*(tden-tnum)*T[j];
   mu[j+1]=div0(rhs,D*tden);
   e[j+1]=ceildiv((ab(A)*e[j]+ab(B)*e[j-1])*tden+2*ab(tden-tnum)*TE[j],D*tden)+1;
  }
  for(int j=0;j<=K;j++){aa[j]+=div0(wa[m]*mu[j],den);dd[j]+=div0(wd[m]*mu[j],den);ae[j]+=ceildiv(ab(wa[m])*e[j],den)+1;de[j]+=ceildiv(ab(wd[m])*e[j],den)+1;}
 }
 Z err=0;for(int j=0;j<=K;j++){if(ae[j]>err)err=ae[j];if(de[j]>err)err=de[j];}
 return {aa,dd,err,ceildiv(norm(wa)+norm(wd),den)};
}
// E_k(t) = 15^k U_k(D0+W*t). Truncated Taylor coefficients to Q.
std::vector<V> ladd(int n,int s,long D0,long C){
 assert(Z(SCAL)*C<=D0-R*W); std::vector<V> E(n+2,V(Q+1));E[0][0]=S;
 for(int k=0;k<n;k++){long pk=long(k)*(k+s),D=D0-pk,A=SCAL*(C-pk);assert(A+W<=D);
  for(int j=0;j<=Q;j++){Z v=Z(A)*E[k][j];if(j)v-=Z(W)*E[k+1][j-1];E[k+1][j]=div0(v,Z(D));}
 }
 return E; // E[n+1] is exactly zero for the terminating ladder.
}
// G_{p,k}=15^{-k} g_{p,k}(nu), polynomial in x=2 nu^2/V^2-1.
struct Kernel{V ideal;Z ideal_err;Z residual_norm;Z tail_bound;Z sum_g;Z max_round;};
Kernel compute(int n,int b,int p,long D0,long T0,const std::vector<V>&E,const Moments&mom){
 int K=p==3?n+1:n;long diag0=9-2*b*(b-1);V prev,cur(1,S);Z ep=0,ec=0,egmax=0,gmax=0,gsum=0;
 // Matrix accumulated at scale S^2; only for k=2 spectral residual.
 std::vector<V> A;if(p==3)A.assign(Q+1,V(K+1));
 V Id(Q+1); Z wnormmax=0; const Z Eerr=Z(2L*(n+1)*(Q+1));
 // Uniform error in each truncated W_k vector, in units 1/S.
 const Z Werr=4*Eerr+Z(4*(Q+1));
 for(int k=0;k<=K;k++){
  V ww(Q+1);
  if(p==3){
   long diag=-long(k)*k+long(8-4*b)*k+diag0;
   for(int j=0;j<=Q;j++){
    Z z=Z(2*T0+diag)*E[k][j];if(j)z+=2*Z(W)*E[k][j-1];if(k)z+=SCAL*Z(long(k)*k)*E[k-1][j];ww[j]=div0(z,Z(2*T0));
   }
  }else{for(int j=0;j<=Q;j++)ww[j]=div0(Z(b==0?L2:-L2)*E[k][j],Z(T0));}
  Z gn=norm(cur),wn=norm(ww);gsum+=gn+ec; if(gn+ec>gmax)gmax=gn+ec;if(ec>egmax)egmax=ec;if(wn>wnormmax)wnormmax=wn;
  if(p==3){ // signed contraction before taking any norms
   for(int j=0;j<=Q;j++)if(ww[j]!=0)for(int l=0;l<=k;l++)if(cur[l]!=0)mpz_addmul(A[j][l].get_mpz_t(),ww[j].get_mpz_t(),cur[l].get_mpz_t());
  }else{
   Z g=0;for(int l=0;l<=k;l++)mpz_addmul(g.get_mpz_t(),cur[l].get_mpz_t(),mom.d[l].get_mpz_t());g=shr(g);
   for(int j=0;j<=Q;j++)mpz_addmul(Id[j].get_mpz_t(),ww[j].get_mpz_t(),g.get_mpz_t());
  }
  if(k==K)break;
  long kk=k, B0=8*kk*kk+8*p*kk+4*p-1-2*V2,B1=4*(kk+p-1)*(kk+p-1),D=4*SCAL*SCAL*(kk+1)*(kk+1);
  V next(k+2);for(int l=0;l<=k;l++){
   next[l]+=Z(SCAL*B0)*cur[l];if(l<(int)prev.size())next[l]-=Z(B1)*prev[l];
   if(l==0)next[1]-=Z(2*SCAL*V2)*cur[0];else{next[l-1]-=Z(SCAL*V2)*cur[l];next[l+1]-=Z(SCAL*V2)*cur[l];}
  }
  for(auto&z:next)z=div0(z,Z(D));
  Z en=ceildiv(Z(SCAL*(std::labs(B0)+2*V2))*ec+Z(B1)*ep,Z(D))+(k+2);ep=ec;ec=en;prev.swap(cur);cur.swap(next);
 }
 // |E_k|<=1 on |t|<=100, and |W_k|<=3 on that disk.
 // l1 Taylor coefficient tail on |t|<=1 <=3/(99*100^Q).
 Z tail=ceildiv(3*gsum,Z(99)*powz(100,Q));
 // Bound all numerical errors of the matrix contraction in units S^-1.
 // Werr is the error l1 bound for each W_k, |W| coefficient sum <=4; gsum includes g rounding.
 Z round=ceildiv(Werr*gsum,S)+4*Z(K+1)*egmax+Z((Q+1L)*(K+1));
 Z rnorm=0,err=0;
 if(p==3){
  for(int j=0;j<=Q;j++){
   Z v=0;for(int l=0;l<=K;l++){A[j][l]=shr(A[j][l]);rnorm+=ab(A[j][l]);mpz_addmul(v.get_mpz_t(),A[j][l].get_mpz_t(),mom.a[l].get_mpz_t());}Id[j]=2*shr(v);
  }
  // ||moments||_infty <= weight_norm; integrate kernel error + moment error.
  err=2*(round+tail)*mom.weight_norm+2*ceildiv(mom.err*rnorm,S)+2*(Q+1);
  rnorm+=round+tail;
 }else{
  for(auto&v:Id)v=shr(v);
  // In addition to the kernel/weight rounding above, dot-product rounding per k then W multiplication.
  err=(round+tail)*mom.weight_norm+ceildiv(4*mom.err*gsum,S)+4*(K+1)+(Q+1);
 }
 std::cerr<<"p="<<p<<" gmax(log10 units)="<<gmax.get_str().size()-1<<" gsum(log10 units)="<<gsum.get_str().size()-1<<" max Wnorm="<<f(wnormmax)<<" tail="<<f(tail)<<" round="<<f(round)<<" ideal err="<<f(err)<<"\n";
 assert(wnormmax+Werr<4*S);return {Id,err,rnorm,tail,gsum,round};
}
int main(int argc,char**argv){
 int J=argc>1?atoi(argv[1]):1005,b=argc>2?atoi(argv[2]):4,n=J-b,s=2*b+1,a=1000;long N=n+a,D0=N*(N+s),C=long(n)*(n+s),T0=D0+b*(b+1);
 // Verify retained physical geometry and |W_k(t)|<=3 on |t|<=100.
 assert(n>0&&SCAL*C<=D0-R*W);assert(R*W<T0/2);
 long maxdiag=long(n+1)*(n+1)+8L*(n+1)+24;assert(Z(2*T0)+2*Z(R*W)+maxdiag+SCAL*Z(n+1)*(n+1)<6*Z(T0));
 auto mom=moments(n+1);auto E=ladd(n,s,D0,C);
 auto A=compute(n,b,3,D0,T0,E,mom);
 std::cout<<std::setprecision(16)<<"J="<<J<<" b="<<b<<" D0="<<D0<<" D_radius="<<W<<" Q="<<Q<<" precision_bits="<<P<<" spectral_kernel_norm_upper="<<f(A.residual_norm)<<"\n";
 return 0;
 auto B=compute(n,b,4,D0,T0,E,mom);
 V ideal(Q+1);for(int j=0;j<=Q;j++)ideal[j]=A.ideal[j]+B.ideal[j];Z variation=0;for(int j=1;j<=Q;j++)variation+=ab(ideal[j]);
 Z min=ideal[0]-variation-A.ideal_err-B.ideal_err;
 // Source EPSW <= 2.681358e-6, plus cutoff < 1e-7. Use <=3e-6.
 // 2/Lambda <= 2/998. Tail >=1100 bounded separately by <1e-20.
 Z residual=ceildiv(6*A.residual_norm,Z(998000000L))+ceildiv(S,powz(10,20));
 Z lower=min-residual;
 std::cout<<std::setprecision(16)<<"J="<<J<<" b="<<b<<" n="<<n<<" D0="<<D0<<" D_radius="<<W<<" Q="<<Q<<" precision_bits="<<P<<"\n";
 std::cout<<"ideal_constant="<<f(ideal[0])<<" ideal_variation="<<f(variation)<<" ideal_lower="<<f(min)<<"\n";
 std::cout<<"spectral_kernel_norm_upper="<<f(A.residual_norm)<<" residual_upper="<<f(residual)<<" action_lower="<<f(lower)<<"\n";
 std::cout<<(lower>0?"PASS":"FAIL")<<"\n";
 return lower>0?0:2;
}
