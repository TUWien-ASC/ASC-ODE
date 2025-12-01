#include <iostream>
#include "autodiff.hpp"

using namespace ASC_ode;

// f(x) = x^2 + 3 sin(x)
template <size_t N>
AutoDiff<N,double> f(const AutoDiff<N,double>& x)
{
    return pow(x, 2) + 3.0 * sin(x);
}

int main()
{
    AutoDiff<1,double> x(1.0);
    x.deriv()[0] = 1.0;  // seed derivative

    auto y = f(x);

    std::cout << "f(1) = " << y.value() << "\n";
    std::cout << "f'(1) (AutoDiff) = " << y.deriv()[0] << "\n";

    auto fnum = [](double z){ return z*z + 3.0*std::sin(z); };

    double h = 1e-6;
    double num = (fnum(1+h) - fnum(1-h)) / (2*h);

    std::cout << "f'(1) (numeric) = " << num << "\n";

    return 0;
}
