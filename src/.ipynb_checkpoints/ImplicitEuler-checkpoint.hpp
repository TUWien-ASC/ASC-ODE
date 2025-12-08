#pragma once
#include "timestepper.hpp"
#include "function_algebra.hpp"
#include "newton.hpp"

class ImplicitEuler : public TimeStepper
{
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter>        m_tau;
    std::shared_ptr<ConstantFunction> m_yold;

public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0))
    {
        m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
        auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
        m_equ = ynew - m_yold - m_tau * rhs;          // y - yold - tau*f(y) = 0
    }

    void doStep(double tau, VectorView<double> y) override
    {
        m_yold->set(y);
        m_tau->set(tau);
        NewtonSolver(m_equ, y);
    }
};