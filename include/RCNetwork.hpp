// include/RCNetwork.hpp
#pragma once
#include "nonlinfunc.hpp"   // repo's NonlinearFunction interface
#include <cmath>

enum class VoltageType { ZERO=0, CONST=1, SIN=2 };

class RCNetwork : public NonlinearFunction
{
private:
  double R;
  double C;
  double V0;
  double omega;
  VoltageType vtype;

public:
  RCNetwork(double R_, double C_, VoltageType vt = VoltageType::CONST,
            double V0_ = 1.0, double omega_ = 0.0)
    : R(R_), C(C_), V0(V0_), omega(omega_), vtype(vt) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }

  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    double v = x(0);
    double theta = x(1);

    double Vt = 0.0;
    if (vtype == VoltageType::ZERO) Vt = 0.0;
    else if (vtype == VoltageType::CONST) Vt = V0;
    else Vt = V0 * std::sin(theta);

    // dv/dt
    f(0) = (1.0/(R*C)) * (Vt - v);

    // dtheta/dt
    f(1) = omega;
  }

  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,0) = -1.0/(R*C);
    if (vtype == VoltageType::SIN)
      df(0,1) = (1.0/(R*C)) * V0 * std::cos(x(1));
    else
      df(0,1) = 0.0;
    df(1,0) = 0.0;
    df(1,1) = 0.0;
  }
};
