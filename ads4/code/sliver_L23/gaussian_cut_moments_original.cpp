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
int mpfr_add(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_sub(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_mul(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_div(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_mul_si(mpfr_ptr,mpfr_srcptr,long,mpfr_rnd_t);int mpfr_div_ui(mpfr_ptr,mpfr_srcptr,unsigned long,mpfr_rnd_t);int mpfr_sqrt(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_exp(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_erfc(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_const_pi(mpfr_ptr,mpfr_rnd_t);int mpfr_mul_z(mpfr_ptr,mpfr_srcptr,mpz_srcptr,mpfr_rnd_t);int mpfr_get_z(mpz_ptr,mpfr_srcptr,mpfr_rnd_t);
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
I erfc(const I&a){I r;mpfr_erfc(r.l.a,a.u.a,MPFR_RNDD);mpfr_erfc(r.u.a,a.l.a,MPFR_RNDU);return r;}
I pi(){I r;mpfr_const_pi(r.l.a,MPFR_RNDD);mpfr_const_pi(r.u.a,MPFR_RNDU);return r;}
struct Seq{I prev,cur,mean,var;Seq(I m,I v,I c):cur(c),mean(m),var(v){}void step(long k){I n=mean*cur+(var*prev)*k;prev=std::move(cur);cur=std::move(n);}};
Q readq(){Z p,q;if(!(std::cin>>p>>q)||q<=0)throw std::runtime_error("input rational");Q r(p,q);r.canonicalize();return r;}
int main(){try{int K,digits,N;if(!(std::cin>>K>>digits>>N)||K<0||K>15000||digits<20||digits>10000||N<1)throw std::runtime_error("settings");PREC=4L*digits+128; // enough extra bits; output widths are verified, never presumed
Z scale;mpz_ui_pow_ui(scale.get_mpz_t(),10,digits);std::vector<I> out(K+1);const I L=sqrt(I(997002))+I(100),sq2pi=sqrt(pi()*2);auto st=std::chrono::steady_clock::now();for(int i=0;i<N;++i){Q c0=readq(),s0=readq(),x0=readq();if(c0<0||s0<=0)throw std::runtime_error("center/width");I c(c0),sig(s0),sig2=sig*sig,D=I(100)+sig2*2,root=sqrt(D),m0=sq2pi*sig*(I(1)-(erfc((L-c)/root)+erfc((L+c)/root))/2ul);
I mean=c/1000ul,var=sig2/1000000ul;I w(x0*(c0==0?Q(1,2):Q(1)));I v=(sig2*100/D)/1000000ul;
auto make=[&](const I&a){I diff=c-a;I e=sqrt(sig2*2/D)*exp(I(0)-diff*diff/D);return Seq((c*100+sig2*a*2)/D/1000ul,v,e);};
Seq minus=make(I(0)-L),plus=make(L);I prev,cur=m0,coef=sig2/1000ul;for(long k=0;k<=2L*K;++k){if(!(k&1))out[k/2]=out[k/2]+w*cur;if(k==2L*K)break;I next=mean*cur+(var*prev)*k+coef*(minus.cur-plus.cur);prev=std::move(cur);cur=std::move(next);minus.step(k);plus.step(k);}
if(i%100==99||i+1==N)std::cerr<<"centers "<<i+1<<"/"<<N<<" seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-st).count()<<'\n';}
for(int j=0;j<=K;++j){V l,u;Z a,b;mpfr_mul_z(l.a,out[j].l.a,scale.get_mpz_t(),MPFR_RNDD);mpfr_mul_z(u.a,out[j].u.a,scale.get_mpz_t(),MPFR_RNDU);mpfr_get_z(a.get_mpz_t(),l.a,MPFR_RNDD);mpfr_get_z(b.get_mpz_t(),u.a,MPFR_RNDU);if(a>b)throw std::runtime_error("inverted enclosure");std::cout<<j<<' '<<a<<' '<<b<<'\n';}return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
