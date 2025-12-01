#include <iostream>
#include <fstream>
#include <memory>
#include "nonlinfunc.hpp"
#include "timestepper.hpp"

using namespace ASC_ode;

// ========================================================
// Mass–Spring System:  y' = f(y)
// y = [x, v]
// x' = v
// v' = -(k/m) * x
// ========================================================
class MassSpring : public NonlinearFunction
{
private:
    double mass;
    double stiffness;

public:
    MassSpring(double m, double k) : mass(m), stiffness(k) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    void evaluate(VectorView<double> x, VectorView<double> f) const override
    {
        f(0) = x(1);
        f(1) = -stiffness / mass * x(0);
    }

    void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override
    {
        df = 0.0;
        df(0, 1) = 1.0;
        df(1, 0) = -stiffness / mass;
    }
};


int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cout << "Usage: ./test_ode T_relative steps method\n";
        std::cout << "Example: ./test_ode 4 100 RK2\n";
        return 1;
    }

    double tend_relative = atof(argv[1]);
    int steps = atoi(argv[2]);
    std::string algorithm = argv[3];

    double tend = tend_relative * M_PI;
    double tau = tend / steps;

    Vector<> y = {1.0, 0.0};
    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

    // ========================================================
    // Choose timestepper
    // ========================================================
    std::unique_ptr<TimeStepper> stepper;

    if (algorithm == "explicit")
        stepper = std::make_unique<ExplicitEuler>(rhs);
    else if (algorithm == "improved")
        stepper = std::make_unique<ImprovedEuler>(rhs);
    else if (algorithm == "implicit")
        stepper = std::make_unique<ImplicitEuler>(rhs);
    else if (algorithm == "CN")
        stepper = std::make_unique<CrankNicolson>(rhs);
    else if (algorithm == "RK2")
        stepper = std::make_unique<RungeKutta2>(rhs);
    else
    {
        std::cout << "Choose method: explicit / improved / implicit / CN / RK2\n";
        return 1;
    }

    // ========================================================
    // Write output file
    // ========================================================
    std::ofstream outfile("output_test_ode.txt");

    std::cout << 0.0 << "  " << y(0) << " " << y(1) << "\n";
    outfile << 0.0 << "  " << y(0) << " " << y(1) << "\n";

    for (int i = 0; i < steps; i++)
    {
        stepper->DoStep(tau, y);

        double t = (i + 1) * tau;
        std::cout << t << "  " << y(0) << " " << y(1) << "\n";
        outfile << t << "  " << y(0) << " " << y(1) << "\n";
    }

    return 0;
}
