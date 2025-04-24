#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <array>

static constexpr double G=6.67430e-11;
static constexpr double THETA=0.5;
static constexpr double DT=0.01;

struct Body{
    double r,th,ph;
    double vr,vth,vph;
    double fr,fth,fph;
    double m;
    Body(double R,double T,double P,double M):r(R),th(T),ph(P),vr(0),vth(0),vph(0),fr(0),fth(0),fph(0),m(M){}
    void reset(){fr=fth=fph=0;}
    void step(){
        vr+=fr/m*DT;
        vth+=fth/m*DT;
        vph+=fph/m*DT;
        r+=vr*DT;
        th+=vth*DT/(r>1e-9?r:1e-9);
        ph+=vph*DT/(r>1e-9?r:1e-9)/std::sin(th>1e-9?th:1e-9);
    }
};

struct Region{
    double r0,r1,t0,t1,p0,p1;
    bool contains(const Body&b)const{
        return b.r>=r0&&b.r<r1&&b.th>=t0&&b.th<t1&&b.ph>=p0&&b.ph<p1;
    }
    Region child(int i)const{
        double rm=(r0+r1)/2,tm=(t0+t1)/2,pm=(p0+p1)/2;
        return {(i&1)?rm:r0,(i&1)?r1:rm,(i&2)?tm:t0,(i&2)?t1:tm,(i&4)?pm:p0,(i&4)?p1:pm};
    }
};

class Node{
    Region reg;
    std::unique_ptr<Body> bptr;
    double m,cx,cy,cz;
    bool leaf;
    std::array<std::unique_ptr<Node>,8> ch;
    static std::array<double,3> sph2cart(double r,double th,double ph){
        return {r*std::sin(th)*std::cos(ph),r*std::sin(th)*std::sin(ph),r*std::cos(th)};
    }
public:
    explicit Node(const Region&r):reg(r),bptr(nullptr),m(0),cx(0),cy(0),cz(0),leaf(true){}
    void insert(Body*bp){
        if(!reg.contains(*bp))return;
        if(leaf&&bptr==nullptr){
            bptr=std::make_unique<Body>(*bp);
            m=bp->m;
            auto c=sph2cart(bp->r,bp->th,bp->ph);
            cx=c[0];cy=c[1];cz=c[2];
        }else{
            if(leaf){
                subdiv();
                place(bptr.get());
                bptr.reset();
                leaf=false;
            }
            place(bp);
            auto c=sph2cart(bp->r,bp->th,bp->ph);
            m+=bp->m;
            cx=(cx*(m-bp->m)+c[0]*bp->m)/m;
            cy=(cy*(m-bp->m)+c[1]*bp->m)/m;
            cz=(cz*(m-bp->m)+c[2]*bp->m)/m;
        }
    }
    void force(Body*bp){
        if(m==0)return;
        auto cb=sph2cart(bp->r,bp->th,bp->ph);
        if(leaf&&bptr&&std::fabs(cb[0]-cx)<1e-12&&std::fabs(cb[1]-cy)<1e-12&&std::fabs(cb[2]-cz)<1e-12) return;
        double dx=cx-cb[0],dy=cy-cb[1],dz=cz-cb[2];
        double dist=std::sqrt(dx*dx+dy*dy+dz*dz)+1e-9;
        double size=reg.r1-reg.r0;
        if(leaf||size/dist<THETA){
            double F=G*bp->m*m/(dist*dist);
            double ux=dx/dist,uy=dy/dist,uz=dz/dist;
            double ur=(cb[0]*ux+cb[1]*uy+cb[2]*uz)/bp->r;
            double ut=(std::cos(bp->th)*(std::cos(bp->ph)*ux+std::sin(bp->ph)*uy)-std::sin(bp->th)*uz);
            double up=(-std::sin(bp->ph)*ux+std::cos(bp->ph)*uy)/std::sin(bp->th);
            bp->fr+=F*ur;
            bp->fth+=F*ut;
            bp->fph+=F*up;
        }else{
            for(auto&c:ch) if(c) c->force(bp);
        }
    }
private:
    void subdiv(){for(int i=0;i<8;++i)ch[i]=std::make_unique<Node>(reg.child(i));}
    void place(Body*bp){for(auto&i:ch)if(i->reg.contains(*bp)){i->insert(bp);return;}}
};

int main(){
    std::vector<Body> bodies;
    bodies.emplace_back(1.0,1.0,1.0,5e10);
    bodies.emplace_back(1.2,1.1,0.9,5e10);
    bodies.emplace_back(0.9,1.3,1.1,5e10);
    Region root{0,2,0,M_PI,0,2*M_PI};
    for(int s=0;s<100;++s){
        Node tree(root);
        for(auto&b:bodies)tree.insert(&b);
        for(auto&b:bodies)b.reset();
        for(auto&b:bodies)tree.force(&b);
        for(auto&b:bodies)b.step();
    }
    for(auto&b:bodies)std::cout<<b.r<<" "<<b.th<<" "<<b.ph<<"\n";
}
