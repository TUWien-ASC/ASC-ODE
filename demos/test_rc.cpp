#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include "timestepper.hpp"
#include "rccircuit.hpp"

using namespace ASC_ode;

int main(int argc, char* argv[])
{
    if(argc < 4){
        std::cout << "Usage: ./test_rc tend steps method\n";
        return 1;
    }

    double tend = atof(argv[1]);
    int steps = atoi(argv[2]);
    std::string method = argv[3];

    double tau = tend / steps;

    // initial state: UC(0) = 0, t(0) = 0
    Vector<> y = {0.0, 0.0};

    // RC circuit: R = 100 Ohm, C = 1e-6 F
    auto rhs = std::make_shared<RCCircuit>(100.0, 1e-6);

    std::unique_ptr<TimeStepper> stepper;

    if(method == "explicit")
        stepper = std::make_unique<ExplicitEuler>(rhs);
    else if(method == "improved")
        stepper = std::make_unique<ImprovedEuler>(rhs);
    else if(method == "implicit")
        stepper = std::make_unique<ImplicitEuler>(rhs);
    else if(method == "CN")
        stepper = std::make_unique<CrankNicolson>(rhs);
    else {
        std::cout << "Choose: explicit / improved / implicit / CN\n";
        return 1;
    }

    // -----------------------------
    // Output file per method:
    // rc_explicit.txt
    // rc_improved.txt
    // rc_implicit.txt
    // rc_CN.txt
    // -----------------------------
    std::string fname = "rc_" + method + ".txt";
    std::ofstream out(fname);

    out << 0 << " " << y(0) << " " << y(1) << "\n";

    for(int i = 0; i < steps; i++){
        stepper->DoStep(tau, y);
        out << (i+1) * tau << " " << y(0) << " " << y(1) << "\n";
    }

    std::cout << "Saved: " << fname << "\n";
    return 0;
}
