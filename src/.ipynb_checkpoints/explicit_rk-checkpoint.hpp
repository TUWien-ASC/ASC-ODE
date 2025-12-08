// src/explicit_rk.hpp
#pragma once
#include <vector>
#include <functional>
#include <cassert>
#include <cmath>

// An explicit Runge-Kutta integrator for autonomous ODE y' = f(t,y) or f(y).
// f: (double t, const std::vector<double>& y, std::vector<double>& out)
class ExplicitRungeKutta {
public:
    ExplicitRungeKutta(const std::vector<std::vector<double>>& a_,
                       const std::vector<double>& b_,
                       const std::vector<double>& c_)
      : a(a_), b(b_), c(c_), s((int)b_.size()) {
        assert((int)c_.size() == s);
        assert((int)a.size() == s);
    }

    void step(const std::function<void(double,const std::vector<double>&,std::vector<double>&)>& f,
              double t, double tau,
              const std::vector<double>& y, std::vector<double>& out) const
    {
        int n = (int)y.size();
        out.assign(n, 0.0);
        std::vector<std::vector<double>> k(s, std::vector<double>(n, 0.0));
        std::vector<double> ytemp(n,0.0);
        for (int j=0;j<s;++j) {
            for (int i=0;i<n;++i) {
                double sum = 0.0;
                for (int l=0;l<j;++l) sum += a[j][l] * k[l][i];
                ytemp[i] = y[i] + tau * sum;
            }
            f(t + c[j]*tau, ytemp, k[j]);
        }
        for (int i=0;i<n;++i) {
            double sum = 0.0;
            for (int j=0;j<s;++j) sum += b[j] * k[j][i];
            out[i] = y[i] + tau * sum;
        }
    }

    int stages() const { return s; }

private:
    std::vector<std::vector<double>> a;
    std::vector<double> b;
    std::vector<double> c;
    int s;
};
