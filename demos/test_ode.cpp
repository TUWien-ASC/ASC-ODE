#include <iostream>

#include <nonlinfunc.hpp>
#include <ode.hpp>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
  }
  
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};


int main()
{
  double tend = 4*M_PI;
  int steps = 100;
  Vector<> y = { 1, 0 };  // initializer list
  auto rhs = std::make_shared<MassSpring>();
  
  // ExplicitEuler stepper(rhs);
  ImplicitEuler stepper(rhs);

  double tau = tend/steps;
  for (int i = 0; i < steps; i++)
  {
     stepper.DoStep(tau, y);
     std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
  }
}
