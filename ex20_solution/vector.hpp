#pragma once
#include <vector>
#include <cassert>
#include <cmath>

/*
  Minimal vector class wrapper for matrix operations.
  BUT we do NOT define VectorD here — the project uses
  `using VectorD = std::vector<double>` from mass_spring.hpp.

  That is the ONLY VectorD alias in the entire exercise.
*/

template<typename T>
struct Vec {
    std::vector<T> d;

    Vec() {}
    explicit Vec(size_t n): d(n) {}
    void resize(size_t n) { d.resize(n); }
    size_t size() const { return d.size(); }

    T& operator[](size_t i) { return d[i]; }
    const T& operator[](size_t i) const { return d[i]; }

    void zero() { for(auto& x : d) x = T(0); }
    void fill(T v) { std::fill(d.begin(), d.end(), v); }
};
