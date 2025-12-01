#ifndef TIMERSTEPPER_HPP
#define TIMERSTEPPER_HPP

#include <functional>
#include <exception>

#include "Newton.hpp"


namespace ASC_ode
{
  
  class TimeStepper
  { 
  protected:
    std::shared_ptr<NonlinearFunction> m_rhs;
  public:
    TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
    virtual ~TimeStepper() = default;
    virtual void DoStep(double tau, VectorView<double> y) = 0;
  };

  class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };

  class ImprovedEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      Vector<double> y_tilde(y.size());
      y_tilde = y + 0.5 *tau * m_vecf;
      this->m_rhs->evaluate(y_tilde, m_vecf);
      y += tau * m_vecf;
    }
  };

  class ImplicitEuler : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
  public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = ynew - m_yold - m_tau * m_rhs;
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(tau);
      NewtonSolver(m_equ, y);
    }
  };
    class RungeKutta2 : public TimeStepper
  {
  public:
    using TimeStepper::TimeStepper;

    void DoStep(double tau, VectorView<double> y) override
    {
      const size_t n = m_rhs->dimX();

      Vector<double> k1(n);
      Vector<double> k2(n);
      Vector<double> ytmp(n);

      // k1 = f(y_n)
      m_rhs->evaluate(y, k1);

      // ytmp = y_n + (tau/2) * k1
      for (size_t i = 0; i < n; ++i)
        ytmp(i) = y(i) + 0.5 * tau * k1(i);

      // k2 = f(y_n + tau/2 * k1)
      m_rhs->evaluate(ytmp, k2);

      // y_{n+1} = y_n + tau * k2
      for (size_t i = 0; i < n; ++i)
        y(i) += tau * k2(i);
    }
  };


  class RungeKutta4 : public TimeStepper
  {
  public:
    using TimeStepper::TimeStepper;

    void DoStep(double tau, VectorView<double> y) override
    {
      const size_t n = m_rhs->dimX();

      Vector<double> k1(n);
      Vector<double> k2(n);
      Vector<double> k3(n);
      Vector<double> k4(n);
      Vector<double> ytmp(n);

      // k1 = f(y_n)
      m_rhs->evaluate(y, k1);

      // k2 = f(y_n + tau/2 * k1)
      for (size_t i = 0; i < n; ++i)
        ytmp(i) = y(i) + 0.5 * tau * k1(i);
      m_rhs->evaluate(ytmp, k2);

      // k3 = f(y_n + tau/2 * k2)
      for (size_t i = 0; i < n; ++i)
        ytmp(i) = y(i) + 0.5 * tau * k2(i);
      m_rhs->evaluate(ytmp, k3);

      // k4 = f(y_n + tau * k3)
      for (size_t i = 0; i < n; ++i)
        ytmp(i) = y(i) + tau * k3(i);
      m_rhs->evaluate(ytmp, k4);

      // y_{n+1} = y_n + tau/6 * (k1 + 2k2 + 2k3 + k4)
      for (size_t i = 0; i < n; ++i)
      {
        y(i) += (tau / 6.0) *
                (k1(i) + 2.0 * k2(i) + 2.0 * k3(i) + k4(i));
      }
    }
  };

  class CrankNicolson : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
  public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs)
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0))
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      auto f_old = std::make_shared<ComposeFunction>(rhs, m_yold);

      m_equ = ynew - m_yold - m_tau * (f_old + m_rhs); 
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(0.5 * tau);
      NewtonSolver(m_equ, y);
    }

  };
}



#endif
