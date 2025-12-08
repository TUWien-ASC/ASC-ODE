#pragma once
#include "timestepper.hpp"

class ImprovedEuler : public TimeStepper
{
    Vector<> m_k1, m_k2, m_temp;

public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs),
          m_k1(rhs->dimF()), m_k2(rhs->dimF()), m_temp(rhs->dimF())
    {}

    void doStep(double tau, VectorView<double> y) override
    {
        m_rhs->evaluate(y, m_k1);                 // k1 = f(y_n)
        m_temp = y + tau * m_k1;                  // predictor
        m_rhs->evaluate(m_temp, m_k2);            // k2 = f(y_n + tau k1)
        y += tau * (0.5 * m_k1 + 0.5 * m_k2);      // corrector
    }
};