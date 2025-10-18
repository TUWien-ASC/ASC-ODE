#ifndef NONLINFUNC_H
#define NONLINFUNC_H

#include <vector.hpp>
#include <matrix.hpp>


namespace ASC_ode
{
  using namespace nanoblas;

  class NonlinearFunction
  {
  public:
    virtual ~NonlinearFunction() = default;
    virtual size_t dimX() const = 0;
    virtual size_t dimF() const = 0;
    virtual void evaluate (VectorView<double> x, VectorView<double> f) const = 0;
    virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const = 0;
  };


  class IdentityFunction : public NonlinearFunction
  {
    size_t n;
  public:
    IdentityFunction (size_t _n) : n(_n) { } 
    size_t dimX() const override { return n; }
    size_t dimF() const override { return n; }
    void evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = x;
    }

    void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
      df.diag() = 1.0;
    }
  };



  class ConstantFunction : public NonlinearFunction
  {
    Vector<> val;
  public:
    ConstantFunction (VectorView<double> _val) : val(_val) { }
    void Set(VectorView<double> _val) { val = _val; }
    VectorView<double> Get() const { return val; }
    size_t dimX() const override { return val.size(); }
    size_t dimF() const override { return val.size(); }
    void evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = val;
    }
    void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
    }
  };

  
  
  class SumFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa, fb;
    double faca, facb;
  public:
    SumFunction (std::shared_ptr<NonlinearFunction> _fa,
                 std::shared_ptr<NonlinearFunction> _fb,
                 double _faca, double _facb)
      : fa(_fa), fb(_fb), faca(_faca), facb(_facb) { }

    size_t dimX() const override { return fa->dimX(); }
    size_t dimF() const override { return fa->dimF(); }
    void evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      fa->evaluate(x, f);
      f *= faca;
      Vector<> tmp(dimF());
      fb->evaluate(x, tmp);
      f += facb*tmp;
    }
    void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      fa->evaluateDeriv(x, df);
      df *= faca;
      Matrix<double> tmp(dimF(), dimX());
      fb->evaluateDeriv(x, tmp);
      df += facb*tmp;
    }
  };


  inline auto operator- (std::shared_ptr<NonlinearFunction> fa, std::shared_ptr<NonlinearFunction> fb)
  {
    return std::make_shared<SumFunction>(fa, fb, 1, -1);
  }

  inline auto operator+ (std::shared_ptr<NonlinearFunction> fa, std::shared_ptr<NonlinearFunction> fb)
  {
    return std::make_shared<SumFunction>(fa, fb, 1, 1);
  }

  
  class ScaleFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa;
    double fac;
  public:
    ScaleFunction (std::shared_ptr<NonlinearFunction> _fa,
                   double _fac)
      : fa(_fa), fac(_fac) { }

    size_t dimX() const override { return fa->dimX(); }
    size_t dimF() const override { return fa->dimF(); }
    void evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      fa->evaluate(x, f);
      f *= fac;

    }
    void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      fa->evaluateDeriv(x, df);
      df *= fac;
    }
  };

  inline auto operator* (double a, std::shared_ptr<NonlinearFunction> f)
  {
    return std::make_shared<ScaleFunction>(f, a);
  }




  // fa(fb)
  class ComposeFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa, fb;
  public:
    ComposeFunction (std::shared_ptr<NonlinearFunction> _fa,
                     std::shared_ptr<NonlinearFunction> _fb)
      : fa(_fa), fb(_fb) { }

    size_t dimX() const override { return fb->dimX(); }
    size_t dimF() const override { return fa->dimF(); }
    void evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      Vector<> tmp(fb->dimF());
      fb->evaluate (x, tmp);
      fa->evaluate (tmp, f);
    }
    void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      Vector<> tmp(fb->dimF());
      fb->evaluate (x, tmp);

      Matrix<double> jaca(fa->dimF(), fa->dimX());
      Matrix<double> jacb(fb->dimF(), fb->dimX());

      fb->evaluateDeriv(x, jacb);
      fa->evaluateDeriv(tmp, jaca);

      df = jaca*jacb;
    }
  };
  
  
  inline auto Compose (std::shared_ptr<NonlinearFunction> fa, std::shared_ptr<NonlinearFunction> fb)
  {
    return make_shared<ComposeFunction> (fa, fb);
  }
  
  class EmbedFunction : public NonlinearFunction
  {
    std::shared_ptr<NonlinearFunction> fa;
    size_t firstx, dimx, firstf, dimf;
    size_t nextx, nextf;
  public:
    EmbedFunction (std::shared_ptr<NonlinearFunction> _fa,
                   size_t _firstx, size_t _dimx,
                   size_t _firstf, size_t _dimf)
      : fa(_fa),
        firstx(_firstx), dimx(_dimx), firstf(_firstf), dimf(_dimf),
        nextx(_firstx+_fa->dimX()), nextf(_firstf+_fa->dimF())
    { }

    size_t dimX() const override { return dimx; }
    size_t dimF() const override { return dimf; }
    void evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      fa->evaluate(x.range(firstx, nextx), f.range(firstf, nextf));
    }
    void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0;
      fa->evaluateDeriv(x.range(firstx, nextx),
                        df.rows(firstf, nextf).cols(firstx, nextx));
    }
  };

  
  class Projector : public NonlinearFunction
  {
    size_t size, first, next;
  public:
    Projector (size_t _size, 
               size_t _first, size_t _next)
      : size(_size), first(_first), next(_next) { }

    size_t dimX() const override { return size; }
    size_t dimF() const override { return size; }
    void evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      f.range(first, next) = x.range(first, next);
    }
    void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
      df.diag().range(first, next) = 1;
    }
  };

  
}

#endif
