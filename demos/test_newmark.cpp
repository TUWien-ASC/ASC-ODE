#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;

class RHS : public NonlinearFunction
{
  size_t dimX() const override { return 1; }
  size_t dimF() const override { return 1; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = -x(0);
  }
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df(0,0) = -1;
  }
};


int main()
{
  double tend = 2*M_PI;
  int steps = 100;
  Vector<> x(1);   x(0) = 1.;
  Vector<> dx(1);  dx(0) = 0;
  auto rhs = std::make_shared<RHS>();
  auto mass = std::make_shared<IdentityFunction>(1);
  solveODE_Newmark(tend, steps, x, dx, rhs, mass,
                   [](double t, VectorView<double> x) { std::cout << "t = " << t << ", x = " << x(0) << std::endl; }
                   );
}
