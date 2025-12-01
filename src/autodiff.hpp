#ifndef AUTODIFF_HPP
#define AUTODIFF_HPP

#include <array>
#include <cmath>
#include <iostream>

namespace ASC_ode
{

// =====================================================
// Variable class (template for independent variable xi)
// =====================================================
template <size_t I, typename T = double>
class Variable
{
private:
    T m_val;

public:
    Variable(T v) : m_val(v) {}
    T value() const { return m_val; }
};

// fallback for normal double
template <typename T = double>
T derivative(T, size_t) { return T(0); }


// =====================================================
// AutoDiff<N,T> class
// =====================================================
template <size_t N, typename T = double>
class AutoDiff
{
private:
    T m_val;
    std::array<T, N> m_deriv;

public:
    AutoDiff() : m_val(0.0), m_deriv{} {}

    AutoDiff(T v) : m_val(v), m_deriv{}
    {
        for (size_t i = 0; i < N; i++)
            m_deriv[i] = derivative(v, i);
    }

    template <size_t I>
    AutoDiff(Variable<I, T> var) : m_val(var.value()), m_deriv{}
    {
        m_deriv[I] = 1.0;
    }

    T value() const { return m_val; }
    std::array<T, N> &deriv() { return m_deriv; }
    const std::array<T, N> &deriv() const { return m_deriv; }
};


// =====================================================
// derivative() for AutoDiff type
// =====================================================
template <size_t N, typename T = double>
auto derivative(const AutoDiff<N, T> &v, size_t index)
{
    return v.deriv()[index];
}


// =====================================================
// Printing operator
// =====================================================
template <size_t N, typename T>
std::ostream &operator<<(std::ostream &os, const AutoDiff<N, T> &ad)
{
    os << "Value: " << ad.value() << ", Deriv: [";
    for (size_t i = 0; i < N; i++)
    {
        os << ad.deriv()[i];
        if (i < N - 1) os << ", ";
    }
    os << "]";
    return os;
}


// =====================================================
// Addition
// =====================================================
template <size_t N, typename T = double>
AutoDiff<N, T> operator+(const AutoDiff<N, T> &a, const AutoDiff<N, T> &b)
{
    AutoDiff<N, T> r(a.value() + b.value());
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = a.deriv()[i] + b.deriv()[i];
    return r;
}

template <size_t N, typename T = double>
AutoDiff<N, T> operator+(T a, const AutoDiff<N, T> &b)
{
    return AutoDiff<N, T>(a) + b;
}

template <size_t N, typename T = double>
AutoDiff<N, T> operator+(const AutoDiff<N, T> &a, T b)
{
    return a + AutoDiff<N, T>(b);
}


// =====================================================
// Subtraction
// =====================================================
template <size_t N, typename T = double>
AutoDiff<N, T> operator-(const AutoDiff<N, T> &a)
{
    AutoDiff<N, T> r(-a.value());
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = -a.deriv()[i];
    return r;
}

template <size_t N, typename T = double>
AutoDiff<N, T> operator-(const AutoDiff<N, T> &a, const AutoDiff<N, T> &b)
{
    AutoDiff<N, T> r(a.value() - b.value());
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = a.deriv()[i] - b.deriv()[i];
    return r;
}


// =====================================================
// Multiplication (3 overloads)
// =====================================================
template <size_t N, typename T = double>
AutoDiff<N, T> operator*(const AutoDiff<N, T> &a, const AutoDiff<N, T> &b)
{
    AutoDiff<N, T> r(a.value() * b.value());
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = a.deriv()[i] * b.value() + a.value() * b.deriv()[i];
    return r;
}

template <size_t N, typename T = double>
AutoDiff<N, T> operator*(T a, const AutoDiff<N, T> &b)
{
    return AutoDiff<N, T>(a) * b;
}

template <size_t N, typename T = double>
AutoDiff<N, T> operator*(const AutoDiff<N, T> &a, T b)
{
    return a * AutoDiff<N, T>(b);
}


// =====================================================
// Division (3 full overloads)
// =====================================================
template <size_t N, typename T = double>
AutoDiff<N, T> operator/(const AutoDiff<N, T> &a, const AutoDiff<N, T> &b)
{
    AutoDiff<N, T> r(a.value() / b.value());
    T b2 = b.value() * b.value();

    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = (a.deriv()[i] * b.value() - a.value() * b.deriv()[i]) / b2;

    return r;
}

template <size_t N, typename T = double>
AutoDiff<N, T> operator/(const AutoDiff<N, T> &a, T c)
{
    AutoDiff<N, T> r(a.value() / c);
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = a.deriv()[i] / c;
    return r;
}

template <size_t N, typename T = double>
AutoDiff<N, T> operator/(T c, const AutoDiff<N, T> &a)
{
    AutoDiff<N, T> r(c / a.value());
    T a2 = a.value() * a.value();

    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = -c * a.deriv()[i] / a2;

    return r;
}


// =====================================================
// Elementary functions
// =====================================================
using std::cos;
using std::exp;
using std::log;
using std::sin;

template <size_t N, typename T = double>
AutoDiff<N, T> sin(const AutoDiff<N, T> &a)
{
    AutoDiff<N, T> r(sin(a.value()));
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = cos(a.value()) * a.deriv()[i];
    return r;
}

template <size_t N, typename T = double>
AutoDiff<N, T> cos(const AutoDiff<N, T> &a)
{
    AutoDiff<N, T> r(cos(a.value()));
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = -sin(a.value()) * a.deriv()[i];
    return r;
}

template <size_t N, typename T = double>
AutoDiff<N, T> exp(const AutoDiff<N, T> &a)
{
    AutoDiff<N, T> r(exp(a.value()));
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = exp(a.value()) * a.deriv()[i];
    return r;
}

template <size_t N, typename T = double>
AutoDiff<N, T> log(const AutoDiff<N, T> &a)
{
    AutoDiff<N, T> r(log(a.value()));
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = a.deriv()[i] / a.value();
    return r;
}

template <size_t N, typename T = double>
AutoDiff<N, T> pow(const AutoDiff<N, T> &a, double p)
{
    AutoDiff<N, T> r(std::pow(a.value(), p));
    for (size_t i = 0; i < N; i++)
        r.deriv()[i] = p * std::pow(a.value(), p - 1) * a.deriv()[i];
    return r;
}

} // namespace ASC_ode

#endif
