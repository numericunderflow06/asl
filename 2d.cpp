#include <iostream>
#include <vector>
#include <cmath>
#include <memory>


double G = 6.67430e-11; 

double theta = 0.5;   

double dt = 0.01;     


struct Body {
    double x, y;     
    double vx, vy;   
    double fx, fy;    
    double mass;

    Body(double x_, double y_, double mass_)
        : x(x_), y(y_), vx(0), vy(0), fx(0), fy(0), mass(mass_) {}

 
    void resetForce() {
        fx = fy = 0.0;
    }


    void update(double dt) {
        vx += fx / mass * dt;
        vy += fy / mass * dt;
        x  += vx * dt;
        y  += vy * dt;
    }
};


struct Quad {
    double xmid, ymid;  
    double length;      

    Quad(double x, double y, double l) : xmid(x), ymid(y), length(l) {}


    bool contains(double x, double y) const {
        return (x >= xmid - length/2) && (x <= xmid + length/2)
            && (y >= ymid - length/2) && (y <= ymid + length/2);
    }

 
    Quad NW() const { return Quad(xmid - length/4, ymid + length/4, length/2); }
    Quad NE() const { return Quad(xmid + length/4, ymid + length/4, length/2); }
    Quad SW() const { return Quad(xmid - length/4, ymid - length/4, length/2); }
    Quad SE() const { return Quad(xmid + length/4, ymid - length/4, length/2); }
};


class BHTree {
private:
    Quad region;
    std::unique_ptr<Body> body;   
    double mass;                   
    double comX, comY;             
    bool isExternal;
    std::unique_ptr<BHTree> NW, NE, SW, SE;

public:
    BHTree(const Quad& quad)
        : region(quad), body(nullptr), mass(0), comX(0), comY(0), isExternal(true) {}


    void insert(Body* b) {
        if (!region.contains(b->x, b->y)) return;

        if (body == nullptr && isExternal) {

            body = std::make_unique<Body>(*b);
            mass = b->mass;
            comX = b->x;
            comY = b->y;
        } else {
            
            if (isExternal) {
               
                subdivide();
        
                placeBody(body.get());
                body.reset();
                isExternal = false;
            }

            placeBody(b);

            mass += b->mass;
            comX = (comX * (mass - b->mass) + b->x * b->mass) / mass;
            comY = (comY * (mass - b->mass) + b->y * b->mass) / mass;
        }
    }


    void updateForce(Body* b) {
        if (mass == 0 || (isExternal && body->x == b->x && body->y == b->y)) return;

        double dx = comX - b->x;
        double dy = comY - b->y;
        double dist = std::sqrt(dx*dx + dy*dy);

        if (isExternal || (region.length / dist) < theta) {

            double F = (G * b->mass * mass) / (dist * dist + 1e-7);
            b->fx += F * dx / dist;
            b->fy += F * dy / dist;
        } else {

            if (NW) NW->updateForce(b);
            if (NE) NE->updateForce(b);
            if (SW) SW->updateForce(b);
            if (SE) SE->updateForce(b);
        }
    }

private:

    void subdivide() {
        NW = std::make_unique<BHTree>(region.NW());
        NE = std::make_unique<BHTree>(region.NE());
        SW = std::make_unique<BHTree>(region.SW());
        SE = std::make_unique<BHTree>(region.SE());
    }


    void placeBody(Body* b) {
        if      (region.NW().contains(b->x, b->y)) NW->insert(b);
        else if (region.NE().contains(b->x, b->y)) NE->insert(b);
        else if (region.SW().contains(b->x, b->y)) SW->insert(b);
        else if (region.SE().contains(b->x, b->y)) SE->insert(b);
    }
};

int main() {

    std::vector<Body> bodies;
    bodies.emplace_back(0.3, 0.5, 5e10);
    bodies.emplace_back(0.7, 0.5, 5e10);
    bodies.emplace_back(0.5, 0.8, 5e10);

    int steps = 100;
    for (int step = 0; step < steps; ++step) {

        Quad rootQuad(0.5, 0.5, 1.0);
        BHTree tree(rootQuad);
        for (auto& b : bodies) {
            tree.insert(&b);
        }

        for (auto& b : bodies) b.resetForce();

        for (auto& b : bodies) tree.updateForce(&b);

        for (auto& b : bodies) b.update(dt);
    }


    for (auto& b : bodies) {
        std::cout << "Body at (" << b.x << ", " << b.y << ")\n";
    }
    return 0;
}
