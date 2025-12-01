#include <iostream>
#include <cmath>
#include "../src/autodiff.hpp"

using namespace ASC_ode;

// Pendulum equation: f(theta) = - (g/L) * sin(theta)
constexpr double g = 9.81;
constexpr double L = 1.0;

int main()
{
    using AD = AutoDiff<1,double>;

    // AutoDiff variable: theta = 0.5 rad
    AD theta(0.5);
    theta.deriv()[0] = 1.0;   // d(theta)/d(theta) = 1

    // Pendulum function
    AD f = -(g/L) * sin(theta);

    // Numerical derivative for comparison
    double h = 1e-6;
    double f_num = -(g/L) * sin(0.5);
    double f_plus = -(g/L) * sin(0.5 + h);
    double df_numeric = (f_plus - f_num) / h;

    std::cout << "Pendulum AutoDiff test\n";
    std::cout << "Value f(theta) = " << f.value() << "\n";
    std::cout << "AutoDiff derivative df/dtheta = " << f.deriv()[0] << "\n";
    std::cout << "Numeric derivative df/dtheta = " << df_numeric << "\n";

    return 0;
}
