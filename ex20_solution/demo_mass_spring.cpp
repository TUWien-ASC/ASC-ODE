#include "mass_spring.hpp"
#include "mss_function.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

// simple integrator (symplectic explicit) just for demo output
int main(){
    // Build a 3-mass chain (like double pendulum with 3 masses)
    int n_masses = 3;
    MassSpringSystem sys(n_masses, 2);
    sys.mass = {1.0, 1.0, 1.0};

    // positions: equilateral triangle for spinning top demo OR chain initial positions
    // We'll build a chain along x-axis with spacing 1.0 and then also create constraints
    double spacing = 1.0;
    VectorD q(sys.dimX()), v(sys.dimX());
    for(int i=0;i<n_masses;++i){
        q[i*2 + 0] = i * spacing; // x
        q[i*2 + 1] = 0.0;         // y
        v[i*2 + 0] = 0.0;
        v[i*2 + 1] = 0.0;
    }

    // add distance constraints between consecutive masses (chain)
    sys.addDistanceConstraint(0,1, spacing);
    sys.addDistanceConstraint(1,2, spacing);

    // set initial velocities: small upward velocity on middle mass to excite
    v[1*2 + 1] = 0.5;

    // Time stepping
    double dt = 0.01;
    double T = 4.0;
    std::ofstream ofs("mass_spring_output.txt");
    ofs << "t ";
    for(int i=0;i<n_masses;++i) ofs << "x" << i << " y" << i << " ";
    for(int i=0;i<n_masses;++i) ofs << "vx" << i << " vy" << i << " ";
    ofs << "\n";

    int steps = int(T/dt);
    for(int s=0;s<=steps;++s){
        double t = s * dt;
        ofs << t << " ";
        for(int i=0;i<n_masses;++i){
            ofs << q[i*2 + 0] << " " << q[i*2 + 1] << " ";
        }
        for(int i=0;i<n_masses;++i){
            ofs << v[i*2 + 0] << " " << v[i*2 + 1] << " ";
        }
        ofs << "\n";

        // compute acceleration via augmented constraint solver
        VectorD qdd, lambda;
        sys.computeAccelerationsConstraint(q, v, qdd, lambda, true, 10.0); // with small stabilization

        // explicit symplectic Euler (simple)
        for(int i=0;i<sys.dimX();++i){
            v[i] += qdd[i] * dt;
            q[i] += v[i] * dt;
        }
    }
    ofs.close();
    std::cout<<"Written mass_spring_output.txt\n";
    return 0;
}
