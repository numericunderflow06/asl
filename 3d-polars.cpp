#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <array>

const double G=6.67430e-11;
const double DT=0.01;
const double THETA=0.5;

struct Body{
    double r,th,ph,vr,vth,vph,fr,fth,fph,m;
    Body(double R,double T,double P,double M):r(R),th(T),ph(P),vr(0),vth(0),vph(0),fr(0),fth(0),fph(0),m(M){}
    void reset(){fr=fth=fph=0;}
    void step(){
        vr  += fr/m*DT;
        vth += fth/m*DT;
        vph += fph/m*DT;
        r   += vr*DT;
        th  += vth*DT/(r>1e-9?r:1e-9);
        ph  += vph*DT/(r>1e-9?r:1e-9)/std::sin(th>1e-9?th:1e-9);
    }
};

struct Region{
    double r0,r1,t0,t1,p0,p1;
    bool contains(const Body&b)const{
        return b.r>=r0 && b.r<r1 && b.th>=t0 && b.th<t1 && b.ph>=p0 && b.ph<p1;
    }
    Region child(int i)const{
        double rm=(r0+r1)*0.5, tm=(t0+t1)*0.5, pm=(p0+p1)*0.5;
        return {(i&1)?rm:r0,(i&1)?r1:rm,(i&2)?tm:t0,(i&2)?t1:tm,(i&4)?pm:p0,(i&4)?p1:pm};
    }
};

inline double sphDist(const Body &b,double R,double T,double P){
    double cosg=std::sin(b.th)*std::sin(T)*std::cos(b.ph-P)+std::cos(b.th)*std::cos(T);
    double d2=b.r*b.r+R*R-2*b.r*R*cosg;
    return std::sqrt(d2+1e-12);
}

class Node{
    Region reg;
    std::array<std::unique_ptr<Node>,8> ch;
    double mass,sr,st,sp;
    bool leaf;
    double cR()const{return sr/mass;}
    double cT()const{return st/mass;}
    double cP()const{return sp/mass;}
public:
    explicit Node(const Region&r):reg(r),mass(0),sr(0),st(0),sp(0),leaf(true){}
    void insert(Body*bp){
        if(!reg.contains(*bp))return;
        sr+=bp->m*bp->r;
        st+=bp->m*bp->th;
        sp+=bp->m*bp->ph;
        mass+=bp->m;
        if(leaf&&mass==bp->m)return;
        if(leaf){sub();leaf=false;}
        for(auto &x:ch) if(x->reg.contains(*bp)){x->insert(bp);break;}
    }
    void force(Body*bp){
        if(mass==0||(leaf&&mass==bp->m))return;
        double R=cR(),T=cT(),P=cP();
        double d=sphDist(*bp,R,T,P);
        double size=reg.r1-reg.r0;
        if(leaf||size/d<THETA){
            double F=G*bp->m*mass/(d*d+1e-9);
            double dR=R-bp->r;
            double dT=T-bp->th;
            double dP=P-bp->ph;
            double s=std::sin(bp->th);
            double compR=dR;
            double compT=bp->r*dT;
            double compP=bp->r*s*dP;
            double norm=std::sqrt(compR*compR+compT*compT+compP*compP)+1e-12;
            bp->fr  += F*compR/norm;
            bp->fth += F*compT/norm;
            bp->fph += F*compP/norm;
        }else{
            for(auto &x:ch) if(x) x->force(bp);
        }
    }
private:
    void sub(){for(int i=0;i<8;++i)ch[i]=std::make_unique<Node>(reg.child(i));}
};

int main(){
    std::vector<Body> bodies;
    bodies.emplace_back(1.0,1.0,1.0,5e10);
    bodies.emplace_back(1.2,1.1,0.9,5e10);
    bodies.emplace_back(0.9,1.3,1.1,5e10);
    Region root{0,2,0,M_PI,0,2*M_PI};
    for(int s=0;s<100;++s){
        Node tree(root);
        for(auto &b:bodies) tree.insert(&b);
        for(auto &b:bodies) b.reset();
        for(auto &b:bodies) tree.force(&b);
        for(auto &b:bodies) b.step();
    }
    for(auto &b:bodies) std::cout<<b.r<<" "<<b.th<<" "<<b.ph<<"\n";
}
