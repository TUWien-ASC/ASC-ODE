#ifndef MASS_SPRING_FUNCTION_HPP
#define MASS_SPRING_FUNCTION_HPP

#include "nonlinfunc.hpp"
#include <vector>
#include <cmath> // For mathematical functions like sin, cos, etc.

class MassSpringFunction : public NonlinearFunction {
public:
    MassSpringFunction(double mass, double spring_constant)
        : m_mass(mass), m_k(spring_constant) {}

    // Evaluate the right-hand side of the ODE system
    std::vector<double> eval(const std::vector<double>& x) const override {
        // x[0] = position (x), x[1] = velocity (v)
        double dx1_dt = x[1];  // dx/dt = v (velocity)
        double dx2_dt = -m_k / m_mass * x[0];  // dv/dt = -k/m * x (acceleration)
        
        return {dx1_dt, dx2_dt};
    }

    // Evaluate the Jacobian matrix of the system
    std::vector<std::vector<double>> jacobian(const std::vector<double>& x) const override {
        // The Jacobian matrix is:
        // [ 0,  1 ]
        // [-k/m, 0]
        
        return {
            {0, 1},
            {-m_k / m_mass, 0}
        };
    }

private:
    double m_mass;  // Mass of the object
    double m_k;     // Spring constant
};

#endif // MASS_SPRING_FUNCTION_HPP
