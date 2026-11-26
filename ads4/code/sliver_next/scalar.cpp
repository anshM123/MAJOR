// Exact rational parameter, directed fixed-point signed-transfer evaluator.
// g++ -O3 -std=c++17 scalar.cpp -lgmpxx -lgmp -o scalar
#include <gmpxx.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
using I=mpz_class;
struct C{I r=0,i=0; C()=default; C(I a,I b=0):r(std::move(a)),i(std::move(b)){} };
C operator+(const C&a,const C&b){return C(a.r+b.r,a.i+b.i);} C operator-(const C&a,const C&b){return C(a.r-b.r,a.i-b.i);} C operator*(const C&a,const C&b){return C(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r);} C operator*(const C&a,const I&b){return C(a.r*b,a.i*b);} I l1(const C&a){return abs(a.r)+abs(a.i);} I sq(const C&a){return a.r*a.r+a.i*a.i;} I isqrt(const I&a){I x;mpz_sqrt(x.get_mpz_t(),a.get_mpz_t());return x;} I ceildiv(const I&a,const I&b){if(a<0||b<=0)throw std::runtime_error("ceil");I x;mpz_cdiv_q(x.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return x;}
C divide(const C&a,const C&b,const I&den){if(den<=0)throw std::runtime_error("zero denom");return C((a.r*b.r+a.i*b.i)/den,(a.i*b.r-a.r*b.i)/den);}
struct Ball {C v; I e=0;}; // error is complex l1, in integer units
struct Step {C a,b,d;I na,nb,nd,ld;Step(C A,C B,C D):a(A),b(B),d(D),na(l1(A)),nb(l1(B)),nd(sq(D)),ld(l1(D)){if(nd==0)throw std::runtime_error("pole");}
 Ball eval(const Ball&u,const Ball&v)const {Ball z;z.v=divide(a*u.v+b*v.v,d,nd);z.e=ceildiv((na*u.e+nb*v.e)*ld,nd)+2;return z;}
};
struct Co {C a,b,be,rn,rd,rn0,rd0,c0,c1;};
// Exact Gaussian rational coefficients with common frequency denominator w.
std::vector<Co> coefficients(int p,int b,int J,const I& nr,const I& ni,const I&w){
 int s=2*b+1,q=b*(b+1);std::vector<Co> co;co.reserve(J);C rn0(1),rd0(1);
 for(int j=0;j<J;++j){I w2=w*w;C be((2*p-1-2*j)*w-2*ni,2*nr);C R((2*j+1)*w+2*ni,-2*nr);C rn=R*R,rd(4*(j+1)*w*((j+1)*w+2*ni),-8*(j+1)*w*nr);
 co.push_back({be*(C(2*(s-1)*w)-be*I(2)),be*(be+C(2*w)),be*(be-C(2*s*w)),rn,rd,rn0,rd0,C((8*q+4*(9-2*b*(b-1)))*w2)-be*(be+C(2*(8-4*b)*w)),be*(be+C(2*(9-4*b)*w))});rn0=rn;rd0=rd;}
 return co;
}
struct Result {I P,Q,Pl,Ql,err;};
Result run(int p,const I&dn,const I&dd,const I&hn,const I&hd,const I&w,const I&S,const std::vector<Co>&co){
 Ball x{C(S),0},xm{C(),0},y{C(),0},ym{C(S),0};C ps,qs;I pe=0,qe=0,ww=w*w;
 for(size_t j=0;j<co.size();++j){const Co&v=co[j];
 if(p==3){Step st((C(8*ww*hd)+v.c0*hn)*v.rd0,(v.c1*hn)*v.rn0,v.rd0*I(8*ww*hd));Ball P=st.eval(x,xm),Q=st.eval(y,ym);ps=ps+P.v;qs=qs+Q.v;pe+=P.e;qe+=Q.e;}
 else{ps=ps+x.v;qs=qs+y.v;pe+=x.e;qe+=y.e;}
 if(j+1==co.size())break;
 Step st((C(4*ww*dn*hd)+v.a*I(hn*dd))*v.rn*v.rd0,v.b*I(hn*dd)*v.rn*v.rn0,(C(4*ww*hd*dd)-v.be*I(hn*dd))*v.rd*v.rd0);
 Ball z=st.eval(x,xm),t=st.eval(y,ym);xm=std::move(x);x=std::move(z);ym=std::move(y);y=std::move(t);
 }
 I pm=isqrt(sq(ps)),qm=isqrt(sq(qs));return {pm+1+pe,ceildiv((qm+1+qe)*dd,dn),pm>pe?I(pm-pe):I(0),qm>qe?I((qm-qe)*dd/dn):I(0),pe>qe?pe:qe};
}
int main(int argc,char**argv){try{int digits=argc>1?std::stoi(argv[1]):245,J=argc>2?std::stoi(argv[2]):1500; if(digits<20||digits>5000||J<1||J>1500)throw std::runtime_error("settings");I S;mpz_ui_pow_ui(S.get_mpz_t(),10,digits);int p,b,pp=-1,bb=-1;I dn,dd,hn,hd,nr,ni,w,pr,pi,pw;std::vector<Co>co;
 // input: p b dnum dden hnum hden nu_real_num nu_imag_num nu_den (all exact integers)
 while(std::cin>>p>>b>>dn>>dd>>hn>>hd>>nr>>ni>>w){if((p!=3&&p!=4)||(b!=0&&b!=4)||dn<=0||dd<=dn||hd<=0||w<=0)throw std::runtime_error("input");if(p!=pp||b!=bb||nr!=pr||ni!=pi||w!=pw){co=coefficients(p,b,J,nr,ni,w);pp=p;bb=b;pr=nr;pi=ni;pw=w;}auto r=run(p,dn,dd,hn,hd,w,S,co);std::cout<<r.P<<' '<<r.Q<<' '<<r.Pl<<' '<<r.Ql<<' '<<r.err<<'\n';}
 if(!std::cin.eof())throw std::runtime_error("parse");return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
