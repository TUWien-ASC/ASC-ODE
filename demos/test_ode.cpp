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
  
  solveODE_IE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { std::cout << t << "  " << y(0) << " " << y(1) << std::endl; });



/*
  auto [a,b] = ComputeABfromC (Gauss3c);
  SolveODE_RK(tend, steps, a, b, Gauss3c, y, rhs, 
              [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });


  cout << "Gauss3c = " << Gauss3c << endl;
  cout << "weights = " << b << endl;
  GaussLegendre (Gauss3c, b);
  cout << "with generic function, c = " << Gauss3c << ", weights = " << b << endl;


  Vector<> Radau(3), RadauWeight(3);
  GaussRadau (Radau, RadauWeight);
  // not sure about weights, comput them via ComputeABfromC
  cout << "Radau = " << Radau << ", weight = " << RadauWeight <<  endl;
 */

 
  /*
  SolveODE_RK(tend, steps, Gauss2a, Gauss2b, Gauss2c, y, rhs, 
              [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });
  */

}
