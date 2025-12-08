#pragma once
#include "matrix.hpp"
#include "vector.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

// Simple dense Gaussian elimination with partial pivoting.

inline bool solveLinearSystem(const Matrix& A_in, const VectorD& b_in, VectorD& x_out) {
    int N = A_in.rows();
    if(A_in.cols() != N) return false;
    if(b_in.size() != N) return false;
    Matrix A = A_in; // copy
    VectorD b = b_in; // copy
    x_out.resize(N);

    // make sure b has contiguous data
    for(int k=0;k<N;++k){
        // pivot
        int piv = k;
        double maxv = std::fabs(A(k,k));
        for(int i=k+1;i<N;++i){
            double v = std::fabs(A(i,k));
            if(v > maxv){ maxv = v; piv = i; }
        }
        if(maxv < 1e-15) return false;
        if(piv != k){
            for(int j=k;j<N;++j) std::swap(A(k,j), A(piv,j));
            std::swap(b[k], b[piv]);
        }
        double akk = A(k,k);
        // normalize row k
        for(int j=k+1;j<N;++j) A(k,j) /= akk;
        b[k] /= akk;
        for(int i=k+1;i<N;++i){
            double factor = A(i,k);
            for(int j=k+1;j<N;++j) A(i,j) -= factor * A(k,j);
            b[i] -= factor * b[k];
        }
    }
    // back substitution
    for(int i=N-1;i>=0;--i){
        double s = b[i];
        for(int j=i+1;j<N;++j) s -= A(i,j) * x_out[j];
        x_out[i] = s; // diagonal normalized
    }
    return true;
}
