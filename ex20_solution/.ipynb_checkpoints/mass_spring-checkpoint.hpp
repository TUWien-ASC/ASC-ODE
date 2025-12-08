#pragma once
// mass_spring.hpp


#include <vector>
#include <cmath>
#include <stdexcept>
#include <cassert>

#ifndef VECTOR_D_DEFINED_BY_REPO
// If your repo already has VectorD, include its header above and undef this block.
using VectorD = std::vector<double>;
#endif

struct DistanceConstraint {
    int i, j;      // mass indices
    double rest;   // rest length
    DistanceConstraint(int ii = 0, int jj = 0, double r = 0.0) : i(ii), j(jj), rest(r) {}
};

class MassSpringSystem {
public:
    // n_masses : number of point masses
    // dim : spatial dimension (2 for the exercise)
    MassSpringSystem(int n_masses = 0, int dim = 2)
      : n(n_masses), dim(dim), mass(n_masses, 1.0)
    {
        if (n < 0) throw std::runtime_error("n_masses must be >= 0");
    }

    // public data (simple)
    std::vector<double> mass;                 // length n
    std::vector<DistanceConstraint> cons;     // distance constraints

    int dimX() const { return n * dim; }
    int numMasses() const { return n; }
    int spaceDim() const { return dim; }

    // add constraint between mass i and j with rest length L
    void addDistanceConstraint(int i, int j, double L) {
        if (i < 0 || j < 0 || i >= n || j >= n) throw std::runtime_error("constraint indices out of range");
        cons.emplace_back(i, j, L);
    }

    // Compute accelerations qdd and Lagrange multipliers lambda for distance constraints.
    // q, v are vectors of length dimX(). qdd and lambda are output (resized inside).
    // stabilization: if true use Baumgarte-like stabilization, stab_k is strength.
    void computeAccelerationsConstraint(const VectorD& q, const VectorD& v,
                                         VectorD& qdd, VectorD& lambda,
                                         bool stabilization = true, double stab_k = 50.0) const
    {
        if ((int)q.size() != dimX() || (int)v.size() != dimX()) {
            throw std::runtime_error("q or v size mismatch");
        }

        // zero right-hand acceleration
        qdd.assign(dimX(), 0.0);

        int m = (int)cons.size();
        if (m == 0) {
            // no constraints -> zero acceleration (no internal forces implemented)
            return;
        }

        // Build J (m x dimX), g (m), Jv (m)
        // We'll compute A = J * M^{-1} * J^T (m x m) and rhs = -J * M^{-1} * f - Jdot*v - alpha*g - beta*Jv
        // Here f = 0 (no springs/forces implemented). We use g = |ri - rj|^2 - rest^2.
        // J rows: for constraint k between i and j: derivative of g w.r.t q is:
        //   dg/dqi = 2*(qi - qj), dg/dqj = -2*(qi - qj)
        // So J_{k, idx} = 2 * (qi - qj) component.

        // Precompute inverse mass diagonal
        VectorD invM(n);
        for (int i = 0; i < n; ++i) {
            double mi = mass[i];
            if (mi <= 0.0) throw std::runtime_error("mass must be > 0");
            invM[i] = 1.0 / mi;
        }

        // Compute g and Jv and store small data for each constraint
        std::vector<double> g(m, 0.0);
        std::vector<double> Jv(m, 0.0);
        // For building A we need J*M^{-1}*J^T. We'll compute A_{kl} by summing contributions explicitly.
        std::vector< std::vector<double> > A(m, std::vector<double>(m, 0.0));

        // store per-constraint vector d = qi - qj (length dim)
        std::vector< std::vector<double> > dvec(m, std::vector<double>(dim, 0.0));

        for (int k = 0; k < m; ++k) {
            int i = cons[k].i;
            int j = cons[k].j;
            double rest = cons[k].rest;

            // d = qi - qj
            for (int dd = 0; dd < dim; ++dd) {
                double qi = q[i*dim + dd];
                double qj = q[j*dim + dd];
                dvec[k][dd] = qi - qj;
            }
            // g = |d|^2 - rest^2
            double dotd = 0.0;
            for (int dd=0; dd<dim; ++dd) dotd += dvec[k][dd] * dvec[k][dd];
            g[k] = dotd - rest * rest;

            // Jv = 2 * d · (vi - vj)
            double dotdv = 0.0;
            for (int dd=0; dd<dim; ++dd) {
                double vi = v[i*dim + dd];
                double vj = v[j*dim + dd];
                dotdv += dvec[k][dd] * (vi - vj);
            }
            Jv[k] = 2.0 * dotdv;
        }

        // Build A matrix (m x m): A_{kl} = (J_k) M^{-1} (J_l)^T
        // Note J_k has components 2*d_k at indices of mass i, and -2*d_k at indices of mass j.
        // So contribution from mass p to A_{kl} = (component at p in J_k)*(component at p in J_l) * invM[p]
        for (int k = 0; k < m; ++k) {
            for (int l = 0; l < m; ++l) {
                double sum = 0.0;
                int ik = cons[k].i, jk = cons[k].j;
                int il = cons[l].i, jl = cons[l].j;

                // contributions from mass ik
                for (int dd=0; dd<dim; ++dd) {
                    double Jk_ik = 2.0 * dvec[k][dd];
                    double Jl_ik = 0.0;
                    if (il == ik) Jl_ik += 2.0 * dvec[l][dd];
                    if (jl == ik) Jl_ik += -2.0 * dvec[l][dd];
                    sum += Jk_ik * Jl_ik * invM[ik];
                }
                // contributions from mass jk
                for (int dd=0; dd<dim; ++dd) {
                    double Jk_jk = -2.0 * dvec[k][dd];
                    double Jl_jk = 0.0;
                    if (il == jk) Jl_jk += 2.0 * dvec[l][dd];
                    if (jl == jk) Jl_jk += -2.0 * dvec[l][dd];
                    sum += Jk_jk * Jl_jk * invM[jk];
                }

                A[k][l] = sum;
            }
        }

        // Build rhs: rhs_k = - Jdot*v - alpha*g - beta*Jv
        // choose simple stabilization coefficients
        double alpha = stabilization ? stab_k : 0.0;
        double beta  = stabilization ? (stab_k * 0.2) : 0.0;

        std::vector<double> rhs(m, 0.0);
        for (int k = 0; k < m; ++k) {
            // Jdot*v is approximated by derivative of Jv; for our g formulation Jv is 2*d·(vi-vj)
            // and Jdot*v approx = 2*(vi-vj)·(vi-vj) + 2*d·(ai - aj) -> we don't have a yet, so we drop second part.
            // We use: rhs = - Jv_dot - alpha*g - beta*Jv. We approximate Jv_dot ~= 2*(vi-vj)·(vi-vj)
            int i = cons[k].i, j = cons[k].j;
            double approx_Jv_dot = 0.0;
            for (int dd=0; dd<dim; ++dd) {
                double dv = v[i*dim+dd] - v[j*dim+dd];
                approx_Jv_dot += dv * dv * 2.0; // scaled
            }
            rhs[k] = - approx_Jv_dot - alpha * g[k] - beta * Jv[k];
        }

        // Solve A * lambda = rhs (m x m). Use naive Gaussian elimination (m is small typically).
        lambda.assign(m, 0.0);
        // If A nearly singular, add small regularization
        double reg = 1e-10;
        for (int i=0;i<m;i++) A[i][i] += reg;

        // Gaussian elimination (in-place)
        std::vector<std::vector<double>> B = A; // copy
        std::vector<double> bb = rhs;
        // forward elimination
        for (int i = 0; i < m; ++i) {
            // pivot
            int piv = i;
            double pivabs = std::fabs(B[i][i]);
            for (int r = i+1; r < m; ++r) {
                if (std::fabs(B[r][i]) > pivabs) { piv = r; pivabs = std::fabs(B[r][i]); }
            }
            if (pivabs < 1e-14) {
                // singular-ish -> set lambda zero (fallback)
                lambda.assign(m, 0.0);
                // compute qdd = - M^{-1} J^T lambda = 0 -> already zero
                return;
            }
            if (piv != i) {
                std::swap(B[i], B[piv]);
                std::swap(bb[i], bb[piv]);
            }
            double diag = B[i][i];
            for (int c = i; c < m; ++c) B[i][c] /= diag;
            bb[i] /= diag;
            for (int r = i+1; r < m; ++r) {
                double fac = B[r][i];
                if (fac == 0.0) continue;
                for (int c = i; c < m; ++c) B[r][c] -= fac * B[i][c];
                bb[r] -= fac * bb[i];
            }
        }
        // back substitution
        for (int i = m-1; i >= 0; --i) {
            double sum = bb[i];
            for (int c = i+1; c < m; ++c) sum -= B[i][c] * lambda[c];
            lambda[i] = sum; // since B[i][i] normalized to 1
        }

        // compute accelerations qdd = M^{-1} * ( - J^T * lambda )
        // initialize qdd zero
        qdd.assign(dimX(), 0.0);
        for (int k = 0; k < m; ++k) {
            int i = cons[k].i;
            int j = cons[k].j;
            double lam = lambda[k];
            // J^T contributes +2*d * lam to mass i, and -2*d * lam to mass j (because J row had +2d for i and -2d for j)
            for (int dd=0; dd<dim; ++dd) {
                double contrib_i = 2.0 * dvec[k][dd] * lam;
                double contrib_j = -2.0 * dvec[k][dd] * lam;
                qdd[i*dim + dd] += contrib_i * invM[i]; // since we multiply by M^{-1} here
                qdd[j*dim + dd] += contrib_j * invM[j];
            }
        }

        // Note: this qdd contains only constraint-induced accelerations (no external forces or springs).
    }

private:
    int n;   // number of masses
    int dim; // spatial dimension
};
