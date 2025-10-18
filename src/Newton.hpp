#ifndef Newton_h
#define Newton_h

#include <lapack_interface.hpp>
#include "nonlinfunc.hpp"

namespace ASC_ode
{
  using namespace std;
  
  void NewtonSolver (std::shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 10,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<double> res(func->dimF());
    Matrix<double> fprime(func->dimF(), func->dimX());

    for (int i = 0; i < maxsteps; i++)
      {
        func->evaluate(x, res);
        double err= norm(res);
        if (err < tol) return;

        func->evaluateDeriv(x, fprime);
        // CalcInverse(fprime);
        // x -= fprime*res;

        LapackLU LU(fprime);
        LU.Solve(res);
        x -= res;
 
        if (callback)
          callback(i, err, x);
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
