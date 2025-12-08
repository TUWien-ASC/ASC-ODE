// demos/demo_autodiff.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>

// ---------- Legendre polynomials P0..Pn ----------
template <typename T>
void LegendrePolynomials(int n, T x, std::vector<T>& P) {
    if (n < 0) { P.clear(); return; }
    P.resize(n+1);
    P[0] = T(1);
    if (n == 0) return;
    P[1] = x;
    for (int k=2;k<=n;++k) {
        P[k] = ((T(2*k-1)*x*P[k-1]) - T(k-1)*P[k-2]) / T(k);
    }
}

// numeric derivative using central difference
double derivative(std::function<double(double)> f, double x, double h=1e-6) {
    return (f(x+h) - f(x-h)) / (2*h);
}

// write Legendre and derivative to CSV
int main_legendre() {
    const int nmax = 5;
    std::ofstream f("legendre.csv");
    f << "x";
    for(int k=0;k<=nmax;++k) f << ",P" << k << ",dP" << k;
    f << "\n";

    const int samples = 201;
    for (int i=0;i<samples;++i) {
        double x = -1.0 + 2.0 * double(i)/(samples-1);
        std::vector<double> P;
        LegendrePolynomials(nmax, x, P);

        f << x;
        for(int k=0;k<=nmax;++k) {
            // numeric derivative
            auto fk = [k](double xx) { std::vector<double> Pk; LegendrePolynomials(5, xx, Pk); return Pk[k]; };
            double dPk = derivative(fk, x);

            f << "," << P[k] << "," << dPk;
        }
        f << "\n";
    }
    f.close();
    std::cout << "Wrote legendre.csv (P0..P5 and derivatives)\n";
    return 0;
}

// ---------- Pendulum test ----------
int main_pendulum() {
    const double g = 9.81;
    const double L = 1.0;
    double alpha = 0.5;       // rad
    double alpha_dot = 0.0;   // rad/s

    // numeric evaluation of derivatives
    auto f0 = [](double alpha, double alpha_dot) { return alpha_dot; };
    auto f1 = [g,L](double alpha, double alpha_dot) { return - (g/L) * sin(alpha); };

    // Jacobian using central difference
    double h = 1e-6;
    double J00 = (f0(alpha+h, alpha_dot) - f0(alpha-h, alpha_dot)) / (2*h);
    double J01 = (f0(alpha, alpha_dot+h) - f0(alpha, alpha_dot-h)) / (2*h);
    double J10 = (f1(alpha+h, alpha_dot) - f1(alpha-h, alpha_dot)) / (2*h);
    double J11 = (f1(alpha, alpha_dot+h) - f1(alpha, alpha_dot-h)) / (2*h);

    std::cout << "Pendulum evaluation at alpha=" << alpha << ", alpha_dot=" << alpha_dot << "\n";
    std::cout << "f0 = " << f0(alpha, alpha_dot) << ", f1 = " << f1(alpha, alpha_dot) << "\n";
    std::cout << "Jacobian (df_i / dx_j):\n";
    std::cout << "[" << J00 << " , " << J01 << "]\n";
    std::cout << "[" << J10 << " , " << J11 << "]\n";
    return 0;
}

int main() {
    main_legendre();
    main_pendulum();
    return 0;
}
