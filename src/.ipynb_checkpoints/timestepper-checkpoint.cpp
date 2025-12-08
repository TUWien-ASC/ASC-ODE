#include "timestepper.hpp"
#include "newton.hpp"

std::vector<double> ExpEuler(const NonlinearFunc& f,
                              const std::vector<double>& x0,
                              double dt) {
    // Explicit Euler method implementation
    std::vector<double> x = x0;
    std::vector<double> fx = f.eval(x);
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += dt * fx[i];
    }
    return x;
}

std::vector<double> ImpEuler(const NonlinearFunc& f,
                              const std::vector<double>& x0,
                              double dt) {
    // Implicit Euler method implementation
    std::vector<double> x = x0;
    
    // Define the system for Newton's method
    auto g = [&](const std::vector<double>& y) -> std::vector<double> {
        std::vector<double> fy = f.eval(y);
        std::vector<double> res(fy.size());
        for (size_t i = 0; i < fy.size(); ++i) {
            res[i] = y[i] - x[i] - dt * fy[i];
        }
        return res;
    };
    
    // Solve using Newton's method
    return NewtonSolver(g, x);
}

std::vector<double> ImprovedEuler(const NonlinearFunc& f,
                                   const std::vector<double>& x0,
                                   double dt) {
    // Improved Euler method implementation
    std::vector<double> x = x0;
    std::vector<double> fx = f.eval(x);
    
    // Take a step with the explicit Euler method
    std::vector<double> x_temp = x;
    for (size_t i = 0; i < x.size(); ++i) {
        x_temp[i] += dt * fx[i];
    }
    
    // Evaluate at the new point
    std::vector<double> fx_temp = f.eval(x_temp);
    
    // Combine the two steps for the improved Euler method
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += dt * 0.5 * (fx[i] + fx_temp[i]);
    }
    
    return x;
}

std::vector<double> CrankNicolson(const NonlinearFunc& f,
                                   const std::vector<double>& x0,
                                   double dt) {
    // Crank-Nicolson method implementation (implicit method)
    std::vector<double> x = x0;
    
    // Define the system for Newton's method
    auto g = [&](const std::vector<double>& y) -> std::vector<double> {
        std::vector<double> fy = f.eval(y);
        std::vector<double> fx = f.eval(x);
        std::vector<double> res(fx.size());
        for (size_t i = 0; i < fx.size(); ++i) {
            res[i] = y[i] - x[i] - 0.5 * dt * (fx[i] + fy[i]);
        }
        return res;
    };
    
    // Solve using Newton's method
    return NewtonSolver(g, x);
}
