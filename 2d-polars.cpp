#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <array>

const double G=6.67430e-11;
const double DT=0.01;
const double THETA=0.5;

struct Body{double r,p,vr,vp,fr,fp,m;Body(double R,double P,double M):r(R),p(P),vr(0),vp(0),fr(0),fp(0),m(M){}void reset(){fr=fp=0;}void step(){vr+=fr/m*DT;vp+=fp/m*DT;r+=vr*DT;p+=vp*DT/(r>1e-9?r:1e-9);} };

struct Region{double r0,r1,p0,p1;bool contains(const Body&b)const{return b.r>=r0&&b.r<r1&&b.p>=p0&&b.p<p1;}Region child(int i)const{double rm=(r0+r1)/2,pm=(p0+p1)/2;return{(i&1)?rm:r0,(i&1)?r1:rm,(i&2)?pm:p0,(i&2)?p1:pm};}};

inline double dist(const Body&b,double R,double P){double d2=b.r*b.r+R*R-2*b.r*R*std::cos(b.p-P);return std::sqrt(d2+1e-12);} 

class Node{Region reg;std::array<std::unique_ptr<Node>,4> ch;double mass,sr,sp;bool leaf;double cR()const{return sr/mass;}double cP()const{return sp/mass;}
public:Node(const Region&r):reg(r),mass(0),sr(0),sp(0),leaf(true){}
 void insert(Body*bp){if(!reg.contains(*bp))return;sr+=bp->m*bp->r;sp+=bp->m*bp->p;mass+=bp->m;if(leaf&&mass==bp->m)return;if(leaf){sub();leaf=false;}for(auto &x:ch)if(x->reg.contains(*bp)){x->insert(bp);break;}}
 void force(Body*bp){if(mass==0||(leaf&&mass==bp->m))return;double R=cR(),P=cP();double d=dist(*bp,R,P);double size=reg.r1-reg.r0;if(leaf||size/d<THETA){double F=G*bp->m*mass/(d*d+1e-9);double dR=R-bp->r;double dP=P-bp->p;double compR=dR;double compP=bp->r*dP;double norm=std::sqrt(compR*compR+compP*compP)+1e-12;bp->fr+=F*compR/norm;bp->fp+=F*compP/norm;}else{for(auto &x:ch)if(x)x->force(bp);}}
private:void sub(){for(int i=0;i<4;++i)ch[i]=std::make_unique<Node>(reg.child(i));}}
;

int main(){std::vector<Body>b; b.emplace_back(1.0,0.2,5e10);b.emplace_back(1.2,1.0,5e10);b.emplace_back(0.8,2.0,5e10);Region root{0,2,0,2*M_PI};for(int s=0;s<100;++s){Node tree(root);for(auto &x:b)tree.insert(&x);for(auto &x:b)x.reset();for(auto &x:b)tree.force(&x);for(auto &x:b)x.step();}for(auto &x:b)std::cout<<x.r<<" "<<x.p<<"\n";}
