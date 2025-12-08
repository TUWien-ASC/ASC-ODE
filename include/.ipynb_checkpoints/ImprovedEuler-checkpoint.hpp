#pragma once
#include "timestepper.hpp"
class ImprovedEuler : public TimeStepper {
    Vector<> k1, k2, tmp;
public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs), k1(rhs->dimF()), k2(rhs->dimF()), tmp(rhs->dimF()) {}
    void doStep(double tau, VectorView<double> y) override {
        m_rhs->evaluate(y, k1);
        tmp = y + tau * k1;
        m_rhs->evaluate(tmp, k2);
        y += tau * (0.5 * k1 + 0.5 * k2);
    }
};
