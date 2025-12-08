// include/timestepper.hpp
#pragma once
#include <memory>
#include <vector>
#include <iostream>
#include <cmath>

// forward declarations from repo's linear algebra
// Replace these with the repo's real Vector / Matrix types.
// Here I show pseudo-types to be adapted to the project's actual types.
#include "nonlinear_function.hpp" // assume this exists in repo (NonlinearFunction)

// The repo uses Vector<> and VectorView<> types; adapt includes if necessary.

class TimeStepper
{
protected:
  std::shared_ptr<NonlinearFunction> m_rhs;
public:
  TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
  virtual ~TimeStepper() = default;
  virtual void doStep(double tau, VectorView<double> y) = 0;
};

/// ExplicitEuler is likely already implemented in the repo; shown for reference
class ExplicitEuler : public TimeStepper
{
  Vector<> m_vecf;
public:
  ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs)
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}

  void doStep(double tau, VectorView<double> y) override
  {
    this->m_rhs->evaluate(y, m_vecf);
    y += tau * m_vecf;
  }
};

/// Improved Euler (Heun / explicit midpoint-like two-stage as defined in exercise)
/// \tilde y = y_n + tau/2 * f(y_n)
/// y_{n+1} = y_n + tau * f( \tilde y )
class ImprovedEuler : public TimeStepper
{
  Vector<> m_vecf;     // f(y)
  Vector<> m_vecft;    // f(tilde y)
  Vector<> m_ytemp;    // temporary vector for tilde y
public:
  ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs)
    : TimeStepper(rhs),
      m_vecf(rhs->dimF()),
      m_vecft(rhs->dimF()),
      m_ytemp(rhs->dimX()) {}

  void doStep(double tau, VectorView<double> y) override
  {
    // evaluate f(y)
    this->m_rhs->evaluate(y, m_vecf);

    // tilde y = y + tau/2 * f(y)
    m_ytemp = y;               // copy current y to temp
    m_ytemp += (tau/2.0) * m_vecf;

    // evaluate f(tilde y)
    this->m_rhs->evaluate(m_ytemp, m_vecft);

    // y_{n+1} = y_n + tau * f(tilde y)
    y += tau * m_vecft;
  }
};

/// CrankNicolson (semi-implicit): we implement a Picard (fixed-point) iteration:
/// y^{k+1} = y_old + 0.5*tau*( f(y_old) + f(y^k) )
/// start from predictor (Improved Euler or Explicit Euler)
class CrankNicolson : public TimeStepper
{
  Vector<> m_fold;     // f(y_old)
  Vector<> m_fk;       // f(y_k)
  Vector<> m_yk;       // current iterate y^k
  int m_max_iter;
  double m_tol;
public:
  CrankNicolson(std::shared_ptr<NonlinearFunction> rhs,
                int max_iter = 20, double tol = 1e-10)
    : TimeStepper(rhs),
      m_fold(rhs->dimF()),
      m_fk(rhs->dimF()),
      m_yk(rhs->dimX()),
      m_max_iter(max_iter),
      m_tol(tol) {}

  void doStep(double tau, VectorView<double> y) override
  {
    // store old
    Vector<> yold = y; // copy (Vector supports copy)

    // evaluate f(yold)
    this->m_rhs->evaluate(yold, m_fold);

    // predictor: one ImprovedEuler step -> starting guess
    // yk = yold + tau * f(yold)  (explicit Euler predictor)
    m_yk = yold;
    m_yk += tau * m_fold;

    // Picard iteration (fixed-point)
    for (int iter = 0; iter < m_max_iter; ++iter)
    {
      // evaluate f(y_k)
      this->m_rhs->evaluate(m_yk, m_fk);

      // compute new iterate
      Vector<> ynew = yold;
      ynew += 0.5 * tau * (m_fold + m_fk);

      // check convergence (norm of difference)
      double diffnorm = 0.0;
      for (size_t i = 0; i < ynew.size(); ++i)
      {
        double d = ynew[i] - m_yk[i];
        diffnorm += d*d;
      }
      diffnorm = std::sqrt(diffnorm);

      m_yk = ynew;

      if (diffnorm < m_tol) break;
    }

    // copy result back to y (view)
    y = m_yk;
  }
};
