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
int mpfr_add(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_sub(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_mul(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_div(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_mul_si(mpfr_ptr,mpfr_srcptr,long,mpfr_rnd_t);int mpfr_div_ui(mpfr_ptr,mpfr_srcptr,unsigned long,mpfr_rnd_t);int mpfr_sqrt(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_exp(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_erfc(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_log(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);int mpfr_const_pi(mpfr_ptr,mpfr_rnd_t);int mpfr_mul_z(mpfr_ptr,mpfr_srcptr,mpz_srcptr,mpfr_rnd_t);int mpfr_get_z(mpz_ptr,mpfr_srcptr,mpfr_rnd_t);
}
#endif
#include <iostream>
#include <fstream>
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
I pi(){I r;mpfr_const_pi(r.l.a,MPFR_RNDD);mpfr_const_pi(r.u.a,MPFR_RNDU);return r;}
I log(const I&a){if(mpfr_cmp_si(a.l.a,0)<=0)throw std::runtime_error("log");I r;mpfr_log(r.l.a,a.l.a,MPFR_RNDD);mpfr_log(r.u.a,a.u.a,MPFR_RNDU);return r;}
Z tenpow(unsigned long n){Z z;mpz_ui_pow_ui(z.get_mpz_t(),10,n);return z;}
Q readq(std::istream&in){Z p,q;if(!(in>>p>>q)||q<=0)throw std::runtime_error("input rational");Q r(p,q);r.canonicalize();return r;}
I pad(const I&v,const I&e){I err;err.u=e.u;return v+err-err;}
void print(const std::string&name,const I&v){Z a,b,S=tenpow(80);V l,u;mpfr_mul_z(l.a,v.l.a,S.get_mpz_t(),MPFR_RNDD);mpfr_mul_z(u.a,v.u.a,S.get_mpz_t(),MPFR_RNDU);mpfr_get_z(a.get_mpz_t(),l.a,MPFR_RNDD);mpfr_get_z(b.get_mpz_t(),u.a,MPFR_RNDU);std::cout<<name<<' '<<a<<' '<<b<<'\n';}
int main(){try{PREC=480;int N;std::ifstream f("quad.txt");f>>N;if(N!=97)throw std::runtime_error("quad count");std::vector<I>x,w;for(int j=0;j<N;++j){x.emplace_back(readq(f));w.emplace_back(readq(f));}
std::vector<I>t,wg;for(int i=0;i<128;i++){I mid(Q(-127+2*i,8));for(int j=0;j<N;j++){I a=mid+x[j]/8ul;t.push_back(a);wg.push_back(w[j]*exp(I(0)-a*a/2ul)/8ul);}}
I pi0=pi(),sq2pi=sqrt(pi0*2),tol(Q(1,tenpow(52)));std::array<I,3> Fs;int cnt=0;std::array<long,3> aa{1,9,49};
int nn;if(!(std::cin>>nn)||nn!=3564)throw std::runtime_error("inventory");for(int i=0;i<nn;i++){Q cq=readq(std::cin),sq=readq(std::cin),aq=readq(std::cin);I c(cq),sig(sq),amp(aq*(cq==0?Q(1,2):Q(1)));for(int j=0;j<3;j++){long a=aa[j];I sum;if(cq<10){if(sq>Q(451,1000))throw std::runtime_error("small c sigma");for(size_t k=0;k<t.size();k++){I nu=c+sig*t[k];sum=sum+wg[k]/(nu*nu*4+I(a));}sum=pad(sum,I(Q(1,tenpow(53))))*sig;}
else{I C0=c*c*4+I(a),C1=c*sig*8,C2=sig*sig*4,fac=C2/C0;I pm,pc=I(1)/C0,gm(1),rho(1);sum=pc;bool pass=false;for(int k=0;k<1000;k++){I pn=(I(0)-C1*pc-C2*pm)/C0;pm=pc;pc=pn;I next=(I(0)-C1*pc-C2*pm)/C0;pm=pc;pc=next;gm=gm*(2*k+1);rho=rho*fac;I err=rho*gm/(unsigned long)a; // degree 2k+1 remainder, before adding term degree 2k+2
if(mpfr_cmp(err.u.a,tol.l.a)<0){sum=pad(sum,err);pass=true;break;}sum=sum+gm*pc;}if(!pass)throw std::runtime_error("Taylor did not pass");sum=sum*sq2pi*sig;}
Fs[j]=Fs[j]+amp*sum;}if(i%500==499)std::cerr<<"pole seeds "<<i+1<<'\n';}
for(int j=0;j<3;j++)print("F"+std::to_string(aa[j]),Fs[j]);std::array<I,3> y;for(int i=0;i<3;i++)y[i]=I(readq(std::cin));I K0,LF9,LF49;
for(int i=0;i<96;i++){I mid(Q(2*i+1,8));for(int j=0;j<N;j++){I nu=mid+x[j]/8ul,z=nu*nu,base=exp(I(0)-z/2ul)*log(I(1)+z)*w[j]/8ul;K0=K0+base;I L=base*(I(1)+z)*(y[0]+z*(y[1]+z*y[2]));LF9=LF9+L/(z*4+I(9));LF49=LF49+L/(z*4+I(49));}}
I tol2(Q(1,tenpow(45)));print("K0",pad(K0,tol2));print("LF9",pad(LF9,tol2));print("LF49",pad(LF49,tol2));print("G0",sqrt(pi0/2ul));print("E0",pi0/2ul*exp(I(Q(1,2)))*erfc(I(1)/sqrt(I(2))));return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
