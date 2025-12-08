#pragma once
#include "timestepper.hpp"
#include "function_algebra.hpp"
#include "newton.hpp"

class CrankNicolson : public TimeStepper
{
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter>        m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
    std::shared_ptr<ConstantFunction> m_fold;   // f(y_n)
    Vector<> m_temp;

public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)),
          m_temp(rhs->dimF())
    {
        m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
        m_fold = std::make_shared<ConstantFunction>(rhs->dimF());
        auto id = std::make_shared<IdentityFunction>(rhs->dimX());

        m_equ = id - m_yold - 0.5 * m_tau * (m_fold + rhs);
    }

    void doStep(double tau, VectorView<double> y) override
    {
        m_rhs->evaluate(y, m_temp);    // compute f(y_n)
        m_fold->set(m_temp);           // store it
        m_yold->set(y);
        m_tau->set(tau);
        NewtonSolver(m_equ, y);
    }
};