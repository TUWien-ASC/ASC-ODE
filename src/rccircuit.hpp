#pragma once
#include "nonlinfunc.hpp"
#include <cmath>

namespace ASC_ode {

class RCCircuit : public NonlinearFunction
{
    double R, C;

public:
    RCCircuit(double R_, double C_) : R(R_), C(C_) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    // x(0) = U_C (capacitor voltage)
    // x(1) = t   (time variable)
    void evaluate(VectorView<double> x, VectorView<double> f) const override {
        double UC = x(0);
        double t  = x(1);
        double U0 = std::cos(100.0 * M_PI * t);

        f(0) = (U0 - UC) / (R*C);   // dUC/dt
        f(1) = 1.0;                 // dt/dt = 1
    }

    void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override {
        df = 0.0;
        df(0,0) = -1.0/(R*C);
        df(0,1) = (-100.0*M_PI * std::sin(100.0 * M_PI * x(1))) / (R*C);
    }
};

}
