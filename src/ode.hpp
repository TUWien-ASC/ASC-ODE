#ifndef ODE_hpp
#define ODE_hpp

#include <functional>
#include <exception>

#include "Newton.hpp"


namespace ASC_ode
{
  
  // implicit Euler method for dy/dt = rhs(y)
  void solveODE_IE(double tend, int steps,
                   VectorView<double> y, std::shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    auto yold = std::make_shared<ConstantFunction>(y);
    auto ynew = std::make_shared<IdentityFunction>(y.size());
    auto equ = ynew-yold - dt * rhs;

    double t = 0;
    for (int i = 0; i < steps; i++)
      {
        NewtonSolver (equ, y);
        yold->set(y);
        t += dt;
        if (callback) callback(t, y);
      }
  }



  


  // implicit Euler method for dy/dt = rhs(y)
  void SolveODE_RK(double tend, int steps,
                   Matrix<> a, Vector<> b, Vector<> c,
                   VectorView<double> y, std::shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    int s = c.size();
    int n = y.size();

    /*
    auto multiple_rhs = make_shared<MultipleFunc>(rhs, s);
    Vector<> my(s*y.Size());
    Vector<> mf(s*y.Size());  
    auto myold = make_shared<ConstantFunction>(my);
    auto mynew = make_shared<IdentityFunction>(s*n);
    auto equ = mynew-myold - dt * Compose(make_shared<MatVecFunc>(a, n), multiple_rhs);
                                          
      
    double t = 0;
    for (int i = 0; i < steps; i++)
      {
        cout << "step " << i << endl;
        for (int j = 0; j < s; j++)
          my.Range(j*n, (j+1)*n) = y;
        myold->Set(my);
        
        NewtonSolver (equ, my);

        multiple_rhs->Evaluate(my, mf);
        for (int j = 0; j < s; j++)
          y += dt * b(j) * mf.Range(j*n, (j+1)*n);
        
        t += dt;
        if (callback) callback(t, y);
      }
    */

    auto multiple_rhs = make_shared<MultipleFunc>(rhs, s);
    Vector<> mk(s*y.size());
    Vector<> my(s*y.size());
    auto myold = std::make_shared<ConstantFunction>(my);
    auto knew = std::make_shared<IdentityFunction>(s*n);
    auto equ = knew - Compose(multiple_rhs, myold+dt*std::make_shared<MatVecFunc>(a, n));
      
    double t = 0;
    for (int i = 0; i < steps; i++)
      {
        std::cout << "step " << i << std::endl;
        for (int j = 0; j < s; j++)
          my.range(j*n, (j+1)*n) = y;
        myold->set(my);

        mk = 0.0;
        NewtonSolver (equ, mk);

        for (int j = 0; j < s; j++)
          y += dt * b(j) * mk.range(j*n, (j+1)*n);
        
        t += dt;
        if (callback) callback(t, y);
      }

  }

  

  
  
  
  
  // Newmark and generalized alpha:
  // https://miaodi.github.io/finite%20element%20method/newmark-generalized/
  
  // Newmark method for  mass*d^2x/dt^2 = rhs
  void solveODE_Newmark(double tend, int steps,
                        VectorView<double> x, VectorView<double> dx,
                        std::shared_ptr<NonlinearFunction> rhs,   
                        std::shared_ptr<NonlinearFunction> mass,  
                        std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    double gamma = 0.5;
    double beta = 0.25;

    Vector<> a(x.size());
    Vector<> v(x.size());

    auto xold = make_shared<ConstantFunction>(x);
    auto vold = make_shared<ConstantFunction>(dx);
    auto aold = make_shared<ConstantFunction>(x);
    rhs->evaluate (xold->get(), aold->get());
    
    auto anew = std::make_shared<IdentityFunction>(a.size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    auto equ = Compose(mass, anew) - Compose(rhs, xnew);

    double t = 0;
    for (int i = 0; i < steps; i++)            
      {
        NewtonSolver (equ, a);
        xnew -> evaluate (a, x);
        vnew -> evaluate (a, v);

        xold->set(x);
        vold->set(v);
        aold->set(a);
        t += dt;
        if (callback) callback(t, x);
      }
    dx = v;
  }




  // Generalized alpha method for M d^2x/dt^2 = rhs
  void SolveODE_Alpha (double tend, int steps, double rhoinf,
                       VectorView<double> x, VectorView<double> dx, VectorView<double> ddx,
                       std::shared_ptr<NonlinearFunction> rhs,   
                       std::shared_ptr<NonlinearFunction> mass,  
                       std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    double alpham = (2*rhoinf-1)/(rhoinf+1);
    double alphaf = rhoinf/(rhoinf+1);
    double gamma = 0.5-alpham+alphaf;
    double beta = 0.25 * (1-alpham+alphaf)*(1-alpham+alphaf);

    Vector<double> a(x.size());
    Vector<double> v(x.size());

    auto xold = make_shared<ConstantFunction>(x);
    auto vold = make_shared<ConstantFunction>(dx);
    auto aold = make_shared<ConstantFunction>(ddx);
    // rhs->Evaluate (xold->Get(), aold->Get()); // solve with M ???
    
    auto anew = std::make_shared<IdentityFunction>(a.size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    // auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - Compose(rhs, (1-alphaf)*xnew+alphaf*xold);
    auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - (1-alphaf)*Compose(rhs,xnew) - alphaf*Compose(rhs, xold);

    double t = 0;
    a = ddx;

    for (int i = 0; i < steps; i++)
      {
        NewtonSolver (equ, a);
        xnew -> evaluate (a, x);
        vnew -> evaluate (a, v);

        xold->set(x);
        vold->set(v);
        aold->set(a);
        t += dt;
        if (callback) callback(t, x);
      }
    dx = v;
    ddx = a;
  }

  

}


#endif
