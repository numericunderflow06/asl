#include <iostream>
#include <vector>
#include <cmath>
#include <memory>

double G = 6.67430e-11;
double theta = 0.5;
double dt = 0.01;

void cartesianToSpherical(double x, double y, double z,
                          double &r, double &theta, double &phi) {
    r = std::sqrt(x*x + y*y + z*z);
    theta = std::atan2(std::sqrt(x*x + y*y), z);
    phi   = std::atan2(y, x);
}

void sphericalToCartesian(double r, double theta, double phi,
                          double &x, double &y, double &z) {
    x = r * std::sin(theta) * std::cos(phi);
    y = r * std::sin(theta) * std::sin(phi);
    z = r * std::cos(theta);
}

struct Body {
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
    double mass;
    double r, theta, phi;

    Body(double x_, double y_, double z_, double mass_)
        : x(x_), y(y_), z(z_), vx(0), vy(0), vz(0),
          fx(0), fy(0), fz(0), mass(mass_) {
        cartesianToSpherical(x, y, z, r, theta, phi);
    }

    void resetForce() {
        fx = fy = fz = 0.0;
    }

    void update(double dt) {
        vx += fx / mass * dt;
        vy += fy / mass * dt;
        vz += fz / mass * dt;
        x  += vx * dt;
        y  += vy * dt;
        z  += vz * dt;
        cartesianToSpherical(x, y, z, r, theta, phi);
    }
};

struct Cube {
    double xmid, ymid, zmid;
    double length;

    Cube(double x, double y, double z, double l)
        : xmid(x), ymid(y), zmid(z), length(l) {}

    bool contains(double x, double y, double z) const {
        return (x >= xmid - length/2 && x <= xmid + length/2) &&
               (y >= ymid - length/2 && y <= ymid + length/2) &&
               (z >= zmid - length/2 && z <= zmid + length/2);
    }

    Cube octant(int i) const {
        double offset = length / 4.0;
        return Cube(
            xmid + offset * ((i&1)? 1:-1),
            ymid + offset * ((i&2)? 1:-1),
            zmid + offset * ((i&4)? 1:-1),
            length / 2.0
        );
    }
};

class BHTree {
private:
    Cube region;
    std::unique_ptr<Body> body;
    double mass;
    double comX, comY, comZ;
    bool isExternal;
    std::array<std::unique_ptr<BHTree>,8> children;

public:
    BHTree(const Cube& cube)
      : region(cube), body(nullptr), mass(0), comX(0), comY(0), comZ(0), isExternal(true) {}

    void insert(Body* b) {
        if (!region.contains(b->x, b->y, b->z)) return;
        if (isExternal && !body) {
            body = std::make_unique<Body>(*b);
            mass = b->mass;
            comX = b->x; comY = b->y; comZ = b->z;
        } else {
            if (isExternal) {
                subdivide();
                placeBody(body.get());
                body.reset();
                isExternal = false;
            }
            placeBody(b);
            mass += b->mass;
            comX = (comX*(mass-b->mass) + b->x*b->mass)/mass;
            comY = (comY*(mass-b->mass) + b->y*b->mass)/mass;
            comZ = (comZ*(mass-b->mass) + b->z*b->mass)/mass;
        }
    }

    void updateForce(Body* b) {
        if (mass == 0 || (isExternal && body && body->x==b->x && body->y==b->y && body->z==b->z)) return;
        double dx = comX - b->x;
        double dy = comY - b->y;
        double dz = comZ - b->z;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz) + 1e-7;
        if (isExternal || (region.length / dist) < theta) {
            double F = G * b->mass * mass / (dist*dist);
            b->fx += F * dx / dist;
            b->fy += F * dy / dist;
            b->fz += F * dz / dist;
        } else {
            for (auto &child : children)
                if (child) child->updateForce(b);
        }
    }

private:
    void subdivide() {
        for (int i=0; i<8; ++i)
            children[i] = std::make_unique<BHTree>(region.octant(i));
    }

    void placeBody(Body* b) {
        for (int i=0; i<8; ++i) {
            if (children[i]->region.contains(b->x, b->y, b->z)) {
                children[i]->insert(b);
                return;
            }
        }
    }
};

int main() {
    std::vector<Body> bodies;
    bodies.emplace_back(0.3, 0.5, 0.2, 5e10);
    bodies.emplace_back(0.7, 0.5, 0.8, 5e10);
    bodies.emplace_back(0.5, 0.8, 0.6, 5e10);
    int steps = 100;
    for (int step=0; step<steps; ++step) {
        Cube root(0.5,0.5,0.5,1.0);
        BHTree tree(root);
        for (auto &b : bodies) tree.insert(&b);
        for (auto &b : bodies) b.resetForce();
        for (auto &b : bodies) tree.updateForce(&b);
        for (auto &b : bodies) b.update(dt);
    }
    for (auto &b : bodies)
        std::cout << "Body at ("<<b.x<<","<<b.y<<","<<b.z<<") r="<<b.r
                  <<" θ="<<b.theta<<" φ="<<b.phi<<"\n";
    return 0;
}
