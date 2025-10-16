#ifndef Newton_h
#define Newton_h

#include <lapack_interface.hpp>
#include "nonlinfunc.hpp"

namespace ASC_ode
{

  void NewtonSolver (std::shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 10,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<double> res(func->DimF());
    Matrix<double> fprime(func->DimF(), func->DimX());

    for (int i = 0; i < maxsteps; i++)
      {
        func->Evaluate(x, res);
        // cout << "|res| = " << L2Norm(res) << endl;
        func->EvaluateDeriv(x, fprime);
        // CalcInverse(fprime);
        // x -= fprime*res;

        auto LU = LapackLU(fprime);
        LU.Solve(res);
        x -= res;
 
        double err= norm(res);
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
