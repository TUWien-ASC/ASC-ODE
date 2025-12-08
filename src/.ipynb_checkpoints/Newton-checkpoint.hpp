#ifndef NEWTON_HPP
#define NEWTON_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "nonlinfunc.hpp"

class NewtonSolver {
public:
    double tol = 1e-8;
    int maxIter = 50;

    void solve(NonlinearFunction& func, std::vector<double>& x) {
        size_t n = x.size();
        std::vector<double> F(n);
        std::vector<std::vector<double>> J(n, std::vector<double>(n));

        for(int iter=0; iter<maxIter; ++iter) {
            func.evaluate(x, F);
            func.evaluateDeriv(x, J);

            // Solve J*dx = -F using simple 2x2 solver
            std::vector<double> dx = linearSolve(J, F);

            for(size_t i=0; i<n; ++i) x[i] -= dx[i];

            double norm = 0.0;
            for(double v : dx) norm += v*v;
            norm = std::sqrt(norm);

            if(norm < tol) {
                std::cout << "Converged in " << iter+1 << " iterations.\n";
                return;
            }
        }
        throw std::runtime_error("Newton did not converge");
    }

private:
    std::vector<double> linearSolve(const std::vector<std::vector<double>>& J, const std::vector<double>& F) {
        // Only works for 2x2 systems
        double a = J[0][0], b = J[0][1];
        double c = J[1][0], d = J[1][1];
        double det = a*d - b*c;
        if(std::abs(det) < 1e-12) throw std::runtime_error("Singular Jacobian");
        return { (d*F[0]-b*F[1])/det, (-c*F[0]+a*F[1])/det };
    }
};

#endif
