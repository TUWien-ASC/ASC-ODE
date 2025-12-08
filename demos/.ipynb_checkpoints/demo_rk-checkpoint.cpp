// demos/demo_rk.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include "src/explicit_rk.hpp"

ExplicitRungeKutta make_explicit_euler() {
    std::vector<std::vector<double>> a = {{0.0}};
    std::vector<double> b = {1.0};
    std::vector<double> c = {0.0};
    return ExplicitRungeKutta(a,b,c);
}
ExplicitRungeKutta make_rk2_midpoint() {
    std::vector<std::vector<double>> a = {{0.0, 0.0},
                                          {0.5, 0.0}};
    std::vector<double> b = {0.0, 1.0};
    std::vector<double> c = {0.0, 0.5};
    return ExplicitRungeKutta(a,b,c);
}
ExplicitRungeKutta make_rk4_classic() {
    std::vector<std::vector<double>> a = {
        {0.0,0.0,0.0,0.0},
        {0.5,0.0,0.0,0.0},
        {0.0,0.5,0.0,0.0},
        {0.0,0.0,1.0,0.0}
    };
    std::vector<double> b = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    std::vector<double> c = {0.0, 0.5, 0.5, 1.0};
    return ExplicitRungeKutta(a,b,c);
}

void pendulum_rhs(double /*t*/, const std::vector<double>& y, std::vector<double>& out) {
    const double g = 9.81, L = 1.0;
    out[0] = y[1];
    out[1] = -(g/L) * std::sin(y[0]);
}

int main() {
    // pendulum compare
    {
        auto euler = make_explicit_euler();
        auto rk2 = make_rk2_midpoint();
        auto rk4 = make_rk4_classic();

        double tend = 10.0;
        int steps = 2000;
        double tau = tend / double(steps);

        std::vector<double> y0 = {0.5, 0.0};

        std::ofstream fe("pendulum_euler.csv");
        std::ofstream f2("pendulum_rk2.csv");
        std::ofstream f4("pendulum_rk4.csv");
        fe << "t,alpha,alpha_dot\n";
        f2 << "t,alpha,alpha_dot\n";
        f4 << "t,alpha,alpha_dot\n";

        std::vector<double> ye = y0, y2 = y0, y4 = y0, ytmp(2);

        for (int i=0;i<=steps;++i) {
            double t = i * tau;
            fe << t << "," << ye[0] << "," << ye[1] << "\n";
            f2 << t << "," << y2[0] << "," << y2[1] << "\n";
            f4 << t << "," << y4[0] << "," << y4[1] << "\n";

            if (i==steps) break;

            euler.step(pendulum_rhs, t, tau, ye, ytmp); ye.swap(ytmp);
            rk2.step(pendulum_rhs, t, tau, y2, ytmp); y2.swap(ytmp);
            rk4.step(pendulum_rhs, t, tau, y4, ytmp); y4.swap(ytmp);
        }
        fe.close(); f2.close(); f4.close();
        std::cout << "Wrote pendulum_{euler,rk2,rk4}.csv\n";
    }

    // linear test for error vs tau
    {
        double lambda = -50.0;
        auto euler = make_explicit_euler();
        auto rk2 = make_rk2_midpoint();
        auto rk4 = make_rk4_classic();

        std::vector<double> taus = {0.02, 0.01, 0.005, 0.0025, 0.00125};
        std::ofstream f("linear_error_vs_tau.csv");
        f << "tau,err_euler,err_rk2,err_rk4\n";

        double tend = 1.0;
        for (double tau : taus) {
            int steps = int(std::ceil(tend / tau));
            double real_tau = tend / double(steps);
            auto exact = [&](double t, double y0){ return y0 * std::exp(lambda * t); };

            double y0 = 1.0;
            std::vector<double> y(1), ytmp(1);

            // Euler
            y = {y0};
            for (int i=0;i<steps;++i) {
                euler.step([&](double t,const std::vector<double>& yy,std::vector<double>& out){
                    out[0] = lambda * yy[0];
                }, i*real_tau, real_tau, y, ytmp);
                y.swap(ytmp);
            }
            double err_euler = std::abs(y[0] - exact(tend,y0));

            // RK2
            y = {y0};
            for (int i=0;i<steps;++i) {
                rk2.step([&](double t,const std::vector<double>& yy,std::vector<double>& out){
                    out[0] = lambda * yy[0];
                }, i*real_tau, real_tau, y, ytmp);
                y.swap(ytmp);
            }
            double err_rk2 = std::abs(y[0] - exact(tend,y0));

            // RK4
            y = {y0};
            for (int i=0;i<steps;++i) {
                rk4.step([&](double t,const std::vector<double>& yy,std::vector<double>& out){
                    out[0] = lambda * yy[0];
                }, i*real_tau, real_tau, y, ytmp);
                y.swap(ytmp);
            }
            double err_rk4 = std::abs(y[0] - exact(tend,y0));

            f << real_tau << "," << err_euler << "," << err_rk2 << "," << err_rk4 << "\n";
        }
        f.close();
        std::cout << "Wrote linear_error_vs_tau.csv\n";
    }

    return 0;
}
