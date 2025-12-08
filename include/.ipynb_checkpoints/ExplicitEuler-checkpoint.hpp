#pragma once
#include "timestepper.hpp"
#include <memory>

// This class implements the standard y_new = y + tau * f(y)
class ExplicitEuler : public TimeStepper {
public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
        : TimeStepper(rhs) {}

    void doStep(double tau, VectorView<double> y) override {
        Vector<> k1; 
        m_rhs->evaluate(y, k1); // Evaluate f(y)
        y += tau * k1;          // Update y
    }
};