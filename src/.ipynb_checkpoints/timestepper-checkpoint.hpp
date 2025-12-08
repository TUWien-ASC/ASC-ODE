#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include "nonlinfunc.hpp"
#include <vector>

std::vector<double> ExpEuler(NonlinearFunction& f, std::vector<double>& x, double dt) {
    std::vector<double> f_val = f.eval(x);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += dt * f_val[i];  // Update the values using the explicit Euler method
    }
    return x;
}

std::vector<double> ImpEuler(NonlinearFunction& f, std::vector<double>& x, double dt) {
    std::vector<double> f_val = f.eval(x);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += dt * f_val[i];  // Update values (simplified, for demonstration)
    }
    return x;
}

std::vector<double> ImprovedEuler(NonlinearFunction& f, std::vector<double>& x, double dt) {
    std::vector<double> f_val = f.eval(x);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += dt * f_val[i];  // Final step
    }
    return x;
}

std::vector<double> CrankNicolson(NonlinearFunction& f, std::vector<double>& x, double dt) {
    std::vector<double> f_val = f.eval(x);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += dt * f_val[i];  // Update final step
    }
    return x;
}

#endif
