#ifndef MYFUNCTION_HPP
#define MYFUNCTION_HPP

#include <vector>
#include "nonlinfunc.hpp"

class MyFunction : public NonlinearFunction {
public:
    double x, y;

    MyFunction(double x0 = 1.0, double y0 = 1.0) : x(x0), y(y0) {}

    // Match the NonlinearFunction interface exactly
    std::vector<double> eval() const {
        return { x*x + y*y - 4, x - y };
    }

    std::vector<std::vector<double>> jacobian() const {
        return {{2*x, 2*y}, {1, -1}};
    }
};

#endif
