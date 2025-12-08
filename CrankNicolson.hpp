#pragma once
#include "timestepper.hpp"
#include "function_algebra.hpp"
#include "newton.hpp"
class CrankNicolson : public TimeStepper {
    std::shared_ptr<NonlinearFunction> eq;
    std::shared_ptr<Parameter> tau;
    std::shared_ptr<ConstantFunction> yold, fold;
    Vector<> ftmp;
public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs), tau(std::make_shared<Parameter>(0.0)), ftmp(rhs->dimF()) {
        yold = std::make_shared<ConstantFunction>(rhs->dimX());
        fold = std::make_shared<ConstantFunction>(rhs->dimF());
        auto id = std::make_shared<IdentityFunction>(rhs->dimX());
        eq = id - yold - 0.5*tau*(fold + rhs);
    }
    void doStep(double tau, VectorView<double> y) override {
        m_rhs->evaluate(y, ftmp); fold->set(ftmp);
        yold->set(y); tau->set(tau);
        NewtonSolver(eq, y);
    }
};
