#pragma once
#include "mass_spring.hpp"
#include "vector.hpp"
#include "matrix.hpp"

// MSS_Function provides the RHS f(x) = [v; qdd] and an exact Jacobian
struct MSS_Function {
    MassSpringSystem sys;

    MSS_Function() = default;
    MSS_Function(const MassSpringSystem& s): sys(s) {}

    // compute RHS: x = [q; v] (size 2*nq)
    void computeRHS(const VectorD& x, VectorD& fx) const;

    // compute exact Jacobian J (size 2*nq x 2*nq)
    void computeJacobian(const VectorD& x, Matrix& J) const;
};
