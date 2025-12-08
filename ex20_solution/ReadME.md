This project implements a **constrained mass–spring system** with the following features:

- Multiple point masses in 2D
- Distance constraints between pairs of masses
- Correct computation of:
  - Positions `q`
  - Velocities `v`
  - Accelerations `q_dd`
  - Lagrange multipliers `λ`
- Augmented constraint solver using a symmetric saddle-point system
- Exact Jacobian `df/dx` for the constrained dynamics
- Demonstration program simulating a 3-mass chain