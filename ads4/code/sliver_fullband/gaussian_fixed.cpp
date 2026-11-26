// Exact-parameter, directed MPFR enclosures for all Gaussian-cutoff moments.
// Compile: g++ -O3 -std=c++17 gaussian_cut_moments.cpp -lgmpxx -lgmp -lmpfr -o gaussian_cut_moments
// (When only the runtime library is installed: replace -lmpfr with -l:libmpfr.so.6.)
#include <gmpxx.h>
#if __has_include(<mpfr.h>)
#include <mpfr.h>
#else
// MPFR 4.x public ABI, for environments that have the runtime but not headers.
extern "C" {
typedef long mpfr_prec_t; typedef long mpfr_exp_t; typedef int mpfr_sign_t;
typedef enum { MPFR_RNDN=0,MPFR_RNDZ,MPFR_RNDU,MPFR_RNDD,MPFR_RNDA,MPFR_RNDF } mpfr_rnd_t;
typedef struct {mpfr_prec_t _mpfr_prec;mpfr_sign_t _mpfr_sign;mpfr_exp_t _mpfr_exp;mp_limb_t *_mpfr_d;} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];typedef __mpfr_struct *mpfr_ptr;typedef const __mpfr_struct *mpfr_srcptr;
void mpfr_init2(mpfr_ptr,mpfr_prec_t); void mpfr_clear(mpfr_ptr); void mpfr_swap(mpfr_ptr,mpfr_ptr);
int mpfr_set(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_set_si(mpfr_ptr,long,mpfr_rnd_t);int mpfr_set_q(mpfr_ptr,mpq_srcptr,mpfr_rnd_t);int mpfr_cmp(mpfr_srcptr,mpfr_srcptr);int mpfr_cmp_si(mpfr_srcptr,long);
int mpfr_add(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_sub(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_mul(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_div(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_mul_2si(mpfr_ptr,mpfr_srcptr,long,mpfr_rnd_t);int mpfr_mul_si(mpfr_ptr,mpfr_srcptr,long,mpfr_rnd_t);int mpfr_div_ui(mpfr_ptr,mpfr_srcptr,unsigned long,mpfr_rnd_t);int mpfr_sqrt(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_exp(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_erfc(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_const_pi(mpfr_ptr,mpfr_rnd_t);int mpfr_mul_z(mpfr_ptr,mpfr_srcptr,mpz_srcptr,mpfr_rnd_t);int mpfr_get_z(mpz_ptr,mpfr_srcptr,mpfr_rnd_t);
}
#endif
#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <chrono>
using Z=mpz_class; using Q=mpq_class; static mpfr_prec_t PREC;
struct V{mpfr_t a;V(){mpfr_init2(a,PREC);mpfr_set_si(a,0,MPFR_RNDN);}V(const V&v):V(){mpfr_set(a,v.a,MPFR_RNDN);}V(V&&v):V(){mpfr_swap(a,v.a);}V&operator=(const V&v){mpfr_set(a,v.a,MPFR_RNDN);return *this;}V&operator=(V&&v){mpfr_swap(a,v.a);return *this;}~V(){mpfr_clear(a);}};
struct I{V l,u;I()=default;I(long n){mpfr_set_si(l.a,n,MPFR_RNDD);mpfr_set_si(u.a,n,MPFR_RNDU);}I(Q q){q.canonicalize();mpfr_set_q(l.a,q.get_mpq_t(),MPFR_RNDD);mpfr_set_q(u.a,q.get_mpq_t(),MPFR_RNDU);} };
I operator+(const I&a,const I&b){I r;mpfr_add(r.l.a,a.l.a,b.l.a,MPFR_RNDD);mpfr_add(r.u.a,a.u.a,b.u.a,MPFR_RNDU);return r;}
I operator-(const I&a,const I&b){I r;mpfr_sub(r.l.a,a.l.a,b.u.a,MPFR_RNDD);mpfr_sub(r.u.a,a.u.a,b.l.a,MPFR_RNDU);return r;}
I operator*(const I&a,const I&b){I r; if(mpfr_cmp_si(a.l.a,0)>=0&&mpfr_cmp_si(b.l.a,0)>=0){mpfr_mul(r.l.a,a.l.a,b.l.a,MPFR_RNDD);mpfr_mul(r.u.a,a.u.a,b.u.a,MPFR_RNDU);return r;}std::array<mpfr_srcptr,2> aa{a.l.a,a.u.a},bb{b.l.a,b.u.a};V x,y;bool first=true;for(auto u:aa)for(auto v:bb){mpfr_mul(x.a,u,v,MPFR_RNDD);mpfr_mul(y.a,u,v,MPFR_RNDU);if(first||mpfr_cmp(x.a,r.l.a)<0)r.l=x;if(first||mpfr_cmp(y.a,r.u.a)>0)r.u=y;first=false;}return r;}
I operator*(const I&a,long b){I r;if(b>=0){mpfr_mul_si(r.l.a,a.l.a,b,MPFR_RNDD);mpfr_mul_si(r.u.a,a.u.a,b,MPFR_RNDU);}else{mpfr_mul_si(r.l.a,a.u.a,b,MPFR_RNDD);mpfr_mul_si(r.u.a,a.l.a,b,MPFR_RNDU);}return r;}
I operator/(const I&a,unsigned long b){if(!b)throw std::runtime_error("division zero");I r;mpfr_div_ui(r.l.a,a.l.a,b,MPFR_RNDD);mpfr_div_ui(r.u.a,a.u.a,b,MPFR_RNDU);return r;}
I operator/(const I&a,const I&b){if(mpfr_cmp_si(b.l.a,0)<=0)throw std::runtime_error("division: denominator not positive");I r;if(mpfr_cmp_si(a.l.a,0)>=0){mpfr_div(r.l.a,a.l.a,b.u.a,MPFR_RNDD);mpfr_div(r.u.a,a.u.a,b.l.a,MPFR_RNDU);return r;}if(mpfr_cmp_si(a.u.a,0)<=0){mpfr_div(r.l.a,a.l.a,b.l.a,MPFR_RNDD);mpfr_div(r.u.a,a.u.a,b.u.a,MPFR_RNDU);return r;}mpfr_div(r.l.a,a.l.a,b.l.a,MPFR_RNDD);mpfr_div(r.u.a,a.u.a,b.l.a,MPFR_RNDU);return r;}
I sqrt(const I&a){if(mpfr_cmp_si(a.l.a,0)<0)throw std::runtime_error("sqrt");I r;mpfr_sqrt(r.l.a,a.l.a,MPFR_RNDD);mpfr_sqrt(r.u.a,a.u.a,MPFR_RNDU);return r;}
I exp(const I&a){I r;mpfr_exp(r.l.a,a.l.a,MPFR_RNDD);mpfr_exp(r.u.a,a.u.a,MPFR_RNDU);return r;}
I erfc(const I&a){I r;
 // Mills: erfc(x)<=exp(-x*x) for x>=0. For a very remote edge
 // retain the rigorous enclosure instead of expensively resolving it.
 // The threshold only affects accuracy, not validity; output widths checked.
 if(mpfr_cmp_si(a.l.a,0)>=0 && mpfr_cmp_si(a.l.a,80)>=0){I t=exp(I(0)-a*a);r.u=t.u;return r;}
 mpfr_erfc(r.l.a,a.u.a,MPFR_RNDD);mpfr_erfc(r.u.a,a.l.a,MPFR_RNDU);return r;}
I scalepow(const I&a,long bits){I r;mpfr_mul_2si(r.l.a,a.l.a,bits,MPFR_RNDD);mpfr_mul_2si(r.u.a,a.u.a,bits,MPFR_RNDU);return r;}
I pi(){I r;mpfr_const_pi(r.l.a,MPFR_RNDD);mpfr_const_pi(r.u.a,MPFR_RNDU);return r;}
// For adjacent exact abscissae x_new<=x_old, integrate exp(-u^2)
// from x_new to x_old using a Taylor series in t in [0,1].
// Cauchy bound on |t|=R: exp(2*x_old*delta*R+delta^2*R^2).
I erfc_increment(const I&old,const I&newx,const I&oldval,const I&tol,const I&twosqrtpi){
 I delta=old-newx;if(mpfr_cmp_si(delta.l.a,0)<0)throw std::runtime_error("nonmonotone center");
 I aa=old*delta*2,bb=delta*delta,front=twosqrtpi*delta*exp(I(0)-old*old);
 int best=100000,lgR=6;I bestr;
 for(int lg=4;lg<=10;lg++){
   long R=1L<<lg;
   I M=twosqrtpi*delta*exp(aa*R+bb*(R*R)-old*old)/(unsigned long)(R-1);
   int m=0;while(mpfr_cmp(scalepow(M,-long(lg)*m).u.a,tol.l.a)>0)m+=8;
   if(m<best){best=m;lgR=lg;bestr=scalepow(M,-long(lg)*m);}
 }
 if(best>20000)throw std::runtime_error("Taylor too long");
 I prev,cur(1),sum(1);
 for(int k=0;k<best;k++){
   I n=(aa*cur-bb*prev*2)/(unsigned long)(k+1);sum=sum+n/(unsigned long)(k+2);prev=std::move(cur);cur=std::move(n);
 }
 I err;mpfr_set(err.u.a,bestr.u.a,MPFR_RNDU);
 I inc=front*sum;return oldval+inc+err-I(0)+I(0)-err; // both-sided error [-bestr,+bestr]
}
Q readq(){Z p,q;if(!(std::cin>>p>>q)||q<=0)throw std::runtime_error("input rational");Q r(p,q);r.canonicalize();return r;}

// All recurrences use nu/1200 and fixed binary scale.  See proof of absolute error.
static long P;static Z S;
Z truncq(const Z&x,const Q&q){Z z=x*q.get_num();mpz_tdiv_q(z.get_mpz_t(),z.get_mpz_t(),q.get_den_mpz_t());return z;}
Z mulfix(const Z&x,const Z&y){Z z=x*y;mpz_tdiv_q_2exp(z.get_mpz_t(),z.get_mpz_t(),P);return z;}
Z mid(const I&a){V l,u;Z x,y;mpfr_mul_2si(l.a,a.l.a,P,MPFR_RNDD);mpfr_mul_2si(u.a,a.u.a,P,MPFR_RNDU);mpfr_get_z(x.get_mpz_t(),l.a,MPFR_RNDD);mpfr_get_z(y.get_mpz_t(),u.a,MPFR_RNDU);if(y<x||y-x>2)throw std::runtime_error("MPFR fixed width >2 units");return (x+y)/2;}
struct FixedSeq{Z prev,cur,mean;Q var;bool skip;FixedSeq(Z m,Q v,Z s):cur(s),mean(m),var(v),skip(s==0){}void step(long k){if(skip)return;Z n=mulfix(cur,mean)+truncq(prev*k,var);prev.swap(cur);cur.swap(n);}};
int main(){try{int K,D,N;if(!(std::cin>>K>>D>>N)||K<1||K>6000||D<100||D>4000||N<1)throw std::runtime_error("settings");Z ten,fac6,fac5;mpz_ui_pow_ui(ten.get_mpz_t(),10,D);mpz_ui_pow_ui(fac6.get_mpz_t(),6,2*K);mpz_ui_pow_ui(fac5.get_mpz_t(),5,2*K);Z target=ten*fac6*Z("100000000000000000000")/fac5;P=mpz_sizeinbase(target.get_mpz_t(),2)+2;S=1;mpz_mul_2exp(S.get_mpz_t(),S.get_mpz_t(),P);PREC=P+128;Z seeds;mpz_ui_pow_ui(seeds.get_mpz_t(),10,800);std::vector<Z> out(K+1);const I L=sqrt(I(997002))+I(100);Q weightMass=0;auto st=std::chrono::steady_clock::now();
for(int i=0;i<N;++i){Q c0=readq(),s0=readq(),x0=readq();Z ml,mu;std::cin>>ml>>mu;if(ml>mu||c0<0||c0>1029||s0<Q(449,1000)||s0>Q(901,1000))throw std::runtime_error("inventory");Q sig2=s0*s0,den=100+2*sig2,mean=c0/1200,var=sig2/1440000,coef=sig2/1200,v=sig2*100/den/1440000,w=x0*(c0==0?Q(1,2):Q(1));weightMass+=abs(w);I c(c0),sg2(sig2),DD(den);
auto make=[&](const I&a){I diff=c-a;Z e=mid(sqrt(sg2*2/DD)*exp(I(0)-diff*diff/DD));return FixedSeq(mid((c*100+sg2*a*2)/DD/1200ul),v,e);};FixedSeq minus=make(I(0)-L),plus=make(L);Z prev,cur=truncq(S,Q((ml+mu)/2,seeds));
for(long k=0;k<=2L*K;++k){if(!(k&1))out[k/2]+=truncq(cur,w);if(k==2L*K)break;Z next=truncq(cur,mean)+truncq(prev*k,var)+truncq(minus.cur-plus.cur,coef);prev.swap(cur);cur.swap(next);minus.step(k);plus.step(k);}
if(i%100==99||i+1==N)std::cerr<<"centers "<<i+1<<"/"<<N<<" seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-st).count()<<'\n';}
if(weightMass>1220000)throw std::runtime_error("mass");Z err=10000000000L, f=ten,g=S;for(int j=0;j<=K;++j){Z a=(out[j]-err)*f,b=(out[j]+err)*f;mpz_fdiv_q(a.get_mpz_t(),a.get_mpz_t(),g.get_mpz_t());mpz_cdiv_q(b.get_mpz_t(),b.get_mpz_t(),g.get_mpz_t());std::cout<<j<<' '<<a<<' '<<b<<'\n';f*=36;g*=25;}
std::cerr<<"PASS fixed scale bits "<<P<<'\n';return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
