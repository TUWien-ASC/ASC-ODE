// mss_function.cpp  
#include "mss_function.hpp"
#include "linear_solve.hpp"
#include <stdexcept>
#include <cmath>

// helper to split x = [q; v]
static void splitState(const VectorD& x, int n_q, VectorD& q, VectorD& v) {
    q.assign(n_q, 0.0);
    v.assign(n_q, 0.0);
    for(int i=0;i<n_q;++i){ q[i]=x[i]; v[i]=x[n_q+i]; }
}

void MSS_Function::computeRHS(const VectorD& x, VectorD& fx) const {
    int n_q = sys.dimX();
    VectorD q, v;
    splitState(x, n_q, q, v);

    VectorD qdd, lambda;
    sys.computeAccelerationsConstraint(q, v, qdd, lambda);

    fx.resize(2*n_q);
    // dx/dt = [ v ; qdd ]
    for(int i=0;i<n_q;++i){ fx[i] = v[i]; fx[n_q + i] = qdd[i]; }
}


// -----------------------------------------------------------------------------
// Exact Jacobian via loop-and-solve:
// df/dx = [ 0    I
//           d(qdd)/dq  d(qdd)/dv ]
//
// We compute qdd by solving augmented system A_aug * x_aug = b_aug where
// x_aug = [ qdd (n_q) ; lambda (nc) ]
// A_aug = [ M     G^T
//           G     0   ]
// b_aug = [ F
//          -Gdotv ]
// -----------------------------------------------------------------------------
void MSS_Function::computeJacobian(const VectorD& x, Matrix& J) const {
    int n_q = sys.dimX();
    int stateN = 2 * n_q;
    J.resize(stateN, stateN); J.zero();

    // top-right block = I
    for(int i=0;i<n_q;++i) J(i, n_q + i) = 1.0;

    // split state
    VectorD q, v;
    splitState(x, n_q, q, v);

    // number of constraints
    int nc = (int)sys.cons.size();

    // Build mass matrix M (n_q x n_q) as diagonal
    Matrix M(n_q, n_q); M.zero();
    for(int mi = 0; mi < sys.numMasses(); ++mi){
        double mval = sys.mass[mi];
        if (mval <= 0.0) throw std::runtime_error("mass must be > 0");
        for(int dd=0; dd<sys.spaceDim(); ++dd){
            int idx = mi*sys.spaceDim() + dd;
            M(idx, idx) = mval;
        }
    }

    // Build G (nc x n_q) and Gdotv (nc) and phi (constraint gap)
    Matrix G(nc, n_q); G.zero();
    VectorD Gdotv(nc, 0.0);
    VectorD phi(nc, 0.0);

    // Precompute dvec per constraint (d = qi - qj)
    std::vector<std::vector<double>> dvec(nc, std::vector<double>(sys.spaceDim(), 0.0));
    for(int r=0;r<nc;++r){
        int i = sys.cons[r].i;
        int j = sys.cons[r].j;
        double rest = sys.cons[r].rest;
        double dotd = 0.0;
        for(int dd=0; dd<sys.spaceDim(); ++dd){
            double di = q[i*sys.spaceDim()+dd];
            double dj = q[j*sys.spaceDim()+dd];
            double dcomp = di - dj;
            dvec[r][dd] = dcomp;
            dotd += dcomp * dcomp;
            // set G at the coordinate positions
            int idx_i = i*sys.spaceDim() + dd;
            int idx_j = j*sys.spaceDim() + dd;
            // Following the constraint formulation used in computeAccelerationsConstraint:
            // G row stores 2*d at i-coordinates and -2*d at j-coordinates
            G(r, idx_i) = 2.0 * dcomp;
            G(r, idx_j) = -2.0 * dcomp;
        }
        phi[r] = dotd - rest*rest;
        // Gdotv = 2 * d · (v_i - v_j)
        double dotdv = 0.0;
        for(int dd=0; dd<sys.spaceDim(); ++dd){
            double vi = v[i*sys.spaceDim()+dd];
            double vj = v[j*sys.spaceDim()+dd];
            dotdv += dvec[r][dd] * (vi - vj);
        }
        Gdotv[r] = 2.0 * dotdv;
    }

    // Build b_aug: top F (no internal forces in simple model -> zero), bottom = -Gdotv
    int Naug = n_q + nc;
    VectorD b(Naug, 0.0);
    // top is F (zero)
    for(int idx=0; idx<n_q; ++idx) b[idx] = 0.0;
    for(int r=0;r<nc;++r) b[n_q + r] = - Gdotv[r];

    // Build augmented matrix A_aug (Naug x Naug)
    Matrix A(Naug, Naug); A.zero();
    // top-left = M
    for(int i=0;i<n_q;++i) for(int j=0;j<n_q;++j) A(i,j) = M(i,j);
    // top-right (i, n_q + r) = G(r,i)
    for(int r=0;r<nc;++r) for(int i=0;i<n_q;++i) A(i, n_q + r) = G(r, i);
    // bottom-left (n_q + r, i) = G(r,i)
    for(int r=0;r<nc;++r) for(int i=0;i<n_q;++i) A(n_q + r, i) = G(r, i);
    // bottom-right zero (already zero)

    // Solve base augmented system A * xaug = b
    VectorD xaug;
    if (!solveLinearSystem(A, b, xaug)) throw std::runtime_error("Base augmented solve failed in Jacobian");

    // Now compute derivative columns
    VectorD rhs(Naug), y(Naug);

    // Precompute xaug components: first n_q are qdd_base, rest are lambda_base (not used directly here)
    // Loop over q-derivatives (d/dq_k)
    for(int k=0;k<n_q;++k){
        // db/dq_k : top part dF/dq_k = 0 (no internal springs here)
        // bottom part d(-Gdotv)/dq_k = - d(Gdotv)/dq_k
        // recall Gdotv[r] = 2 * dvec[r] · (v_i - v_j)
        // derivative wrt a coordinate q_k affects dvec entries:
        // if k corresponds to coordinate dd of mass i -> d(dvec[r][dd])/dq_k = 1 (if that constraint uses i)
        // if k is coordinate dd of mass j -> derivative is -1
        for(int i=0;i<Naug;++i) rhs[i] = 0.0;

        for(int r=0;r<nc;++r){
            const auto &C = sys.cons[r];
            int ii = C.i, jj = C.j;
            for(int dd=0; dd<sys.spaceDim(); ++dd){
                int idx_i = ii*sys.spaceDim() + dd;
                int idx_j = jj*sys.spaceDim() + dd;
                if(k == idx_i){
                    // d(Gdotv)/dq_k contribution: 2 * (v_i - v_j)[dd]
                    rhs[n_q + r] -= 2.0 * ( v[ii*sys.spaceDim()+dd] - v[jj*sys.spaceDim()+dd] );
                } else if(k == idx_j){
                    // derivative of dvec is -1 -> contribution -2*(v_i - v_j) * (-1) => +2*(v_i - v_j)
                    rhs[n_q + r] -= 2.0 * ( - (v[ii*sys.spaceDim()+dd] - v[jj*sys.spaceDim()+dd]) );
                    // simplify: that's +2*(v[ii]-v[jj])? but above sign covers
                }
            }
            // Note: we DID NOT include stabilization derivatives here; if you use Baumgarte terms, add them.
        }

        // Compute (dA/dq_k) * xaug and subtract: rhs_total = db_dqk - dA_dqk * xaug
        // dA/dq_k nonzero only in entries where G depends on q (top-right and bottom-left blocks)
        // For a given constraint r and coordinate c:
        // dG(r,c)/dq_k = derivative of 2*d_component: equals 2 if c==coord of mass i (same dd) and k==that coord,
        // equals -2 if c==coord of mass j and k==that coord, else 0.
        // Implement multiply dA * xaug without building dA.

        VectorD temp(Naug, 0.0); // temp = (dA/dq_k) * xaug
        for(int r=0;r<nc;++r){
            const auto &C = sys.cons[r];
            int ii = C.i, jj = C.j;
            for(int dd=0; dd<sys.spaceDim(); ++dd){
                int idx_i = ii*sys.spaceDim() + dd;
                int idx_j = jj*sys.spaceDim() + dd;

                if(k == idx_i){
                    // dG(r, idx_i) / dq_k = 2
                    // contributes to top-row idx_i: (dA)_{idx_i, n_q + r} = 2 -> multiply by xaug[n_q + r]
                    temp[idx_i] += 2.0 * xaug[n_q + r];
                    // contributes to bottom-row n_q + r: (dA)_{n_q + r, idx_i} = 2 -> multiply by xaug[idx_i]
                    temp[n_q + r] += 2.0 * xaug[idx_i];
                }
                if(k == idx_j){
                    // dG(r, idx_j) / dq_k = -2
                    temp[idx_j] += -2.0 * xaug[n_q + r];
                    temp[n_q + r] += -2.0 * xaug[idx_j];
                }
            }
        }

        // subtract temp
        for(int i=0;i<Naug;++i) rhs[i] -= temp[i];

        // solve A * y = rhs
        if(!solveLinearSystem(A, rhs, y)) throw std::runtime_error("Jacobian column solve failed (dq)");
        // store accelerations derivative dqdd/dq_k in Jacobian bottom-left block (rows n_q.., col k)
        for(int i=0;i<n_q;++i) J(n_q + i, k) = y[i];
    }

    // now derivatives w.r.t. velocities (d/dv_k)
    for(int k=0;k<n_q;++k){
        for(int i=0;i<Naug;++i) rhs[i] = 0.0;

        // db/dv_k : top part zero, bottom part = - d(Gdotv)/dv_k
        // Gdotv[r] = 2 * sum_c (dvec[r,c] * v_c_contribution)
        // derivative wrt v_k is simply 2 * d_component located at column k (i.e. G(r,k))
        for(int r=0;r<nc;++r){
            rhs[n_q + r] -= G(r, k);
        }

        // dA/dv_k = 0 (G does not depend on v), so no dA term
        if(!solveLinearSystem(A, rhs, y)) throw std::runtime_error("Jacobian column solve failed (dv)");
        for(int i=0;i<n_q;++i) J(n_q + i, n_q + k) = y[i];
    }

    // Top-left is zero already; top-right was set above; full Jacobian ready.
}
