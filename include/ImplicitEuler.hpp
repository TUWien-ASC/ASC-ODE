#pragma once
#include "timestepper.hpp"
#include "function_algebra.hpp"
#include "newton.hpp"
class ImplicitEuler : public TimeStepper {
    std::shared_ptr<NonlinearFunction> eq;
    std::shared_ptr<Parameter> tau;
    std::shared_ptr<ConstantFunction> yold;
public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs), tau(std::make_shared<Parameter>(0.0)) {
        yold = std::make_shared<ConstantFunction>(rhs->dimX());
        auto id = std::make_shared<IdentityFunction>(rhs->dimX());
        eq = id - yold - tau * rhs;
    }
    void doStep(double tau, VectorView<double> y) override {
        yold->set(y); tau->set(tau); NewtonSolver(eq, y);
    }
};
