// include/mass_spring.hpp
#pragma once
#include "nonlinear_function.hpp"  // repo's interface

class MassSpring : public NonlinearFunction
{
  double m_mass;
  double m_stiffness;

public:
  MassSpring (double m, double k)
    : m_mass(m), m_stiffness(k) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }

  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -m_stiffness/m_mass * x(0);
  }

  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -m_stiffness/m_mass;
  }
};
