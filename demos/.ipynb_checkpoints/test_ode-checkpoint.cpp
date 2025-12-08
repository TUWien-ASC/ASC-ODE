// demos/test_ode.cpp
#include <iostream>
#include <fstream>
#include <memory>

#include "mass_spring.hpp"   // implement the MassSpring NonlinearFunction (see below)
#include "timestepper.hpp"   // include the file we modified/added

int main()
{
  double m = 1.0, k = 1.0;
  double tend = 4*M_PI;
  int steps = 1000;   // try different time-steps to see behavior
  double tau = tend / double(steps);

  // initial y = [position, velocity]
  Vector<> y0 = { 1.0, 0.0 };

  // create rhs
  std::shared_ptr<NonlinearFunction> rhs = std::make_shared<MassSpring>(m, k);

  // explicit Euler
  {
    ExplicitEuler ee(rhs);
    Vector<> y = y0;
    std::ofstream f("mass_explicit.csv");
    f << "t,x,v\n";
    f << 0.0 << "," << y(0) << "," << y(1) << "\n";
    for (int i = 0; i < steps; ++i)
    {
      ee.doStep(tau, y);
      f << (i+1)*tau << "," << y(0) << "," << y(1) << "\n";
    }
    f.close();
    std::cout << "Wrote mass_explicit.csv\n";
  }

  // improved Euler
  {
    ImprovedEuler ie(rhs);
    Vector<> y = y0;
    std::ofstream f("mass_improved.csv");
    f << "t,x,v\n";
    f << 0.0 << "," << y(0) << "," << y(1) << "\n";
    for (int i = 0; i < steps; ++i)
    {
      ie.doStep(tau, y);
      f << (i+1)*tau << "," << y(0) << "," << y(1) << "\n";
    }
    f.close();
    std::cout << "Wrote mass_improved.csv\n";
  }

  // crank-nicolson (picard)
  {
    CrankNicolson cn(rhs, 50, 1e-12);
    Vector<> y = y0;
    std::ofstream f("mass_crank.csv");
    f << "t,x,v\n";
    f << 0.0 << "," << y(0) << "," << y(1) << "\n";
    for (int i = 0; i < steps; ++i)
    {
      cn.doStep(tau, y);
      f << (i+1)*tau << "," << y(0) << "," << y(1) << "\n";
    }
    f.close();
    std::cout << "Wrote mass_crank.csv\n";
  }

  std::cout << "Done.\n";
  return 0;
}
