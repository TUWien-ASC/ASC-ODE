# Numerical Methods for ODEs — Homework Report  


---

## 1. Introduction
This report summarizes the implementation and results of several numerical methods for solving ordinary differential equations (ODEs). It covers:

- **HW 17.4:** Implicit Euler, Crank–Nicolson, and RC circuit simulation  
- **HW 18:** Forward-mode Automatic Differentiation  
- **HW 19:** Runge–Kutta methods (RK2)

All methods were integrated into a unified `TimeStepper` framework and tested using the provided mass–spring and RC circuits.

---

## 2. HW 17.4 — Implicit Euler, Crank–Nicolson & RC Circuit

### 2.1 Mass–Spring System
We study the classical system:

$$
x' = v, \qquad v' = -\frac{k}{m}x
$$

This describes a harmonic oscillator.  
The following methods were implemented and tested:

- Explicit Euler  
- Improved Euler  
- **Implicit Euler**  
- **Crank–Nicolson**

#### **Plots**
Plots for the solution trajectories and phase-space behavior were generated.
<h4 style="text-align:center;">Time Evolution (Mass–Spring)</h4>
<p align="center">
  <img src="report_images/Time_evolution_standard.png" width="450">
</p>

<h4 style="text-align:center;">Phase Plot (Mass–Spring)</h4>
<p align="center">
  <img src="report_images/Phase_plot_standard.png" width="450">
</p>

#### **Observations**
- **Explicit Euler**: noticeable phase error and amplitude drift.  
- **Improved Euler**: better accuracy but still accumulates drift.  
- **Implicit Euler**: unconditionally stable, slightly overdamped.  
- **Crank–Nicolson**: best accuracy among Euler-type schemes, symmetric and stable.

---

### 2.2 RC Circuit ODE
The charging of a capacitor follows:

$$
C \frac{dU_c}{dt} = \frac{U_s - U_c}{R}
$$

which was written in autonomous form and solved using all four methods.
#### **Plots**

<h4 style="text-align:center;">RC Circuit – Full Time Evolution</h4>
<p align="center">
  <img src="report_images/RC_Time_evolution_tend_0.01_steps_2000.png" width="500">
</p>


<h4 style="text-align:center;">RC Circuit – Zoomed View (Stability)</h4>
<p align="center">
  <img src="report_images/RC_Time_evolution_zoom_tend_0.01_steps_2000.png" width="500">
</p>

#### **Observations**
- Explicit Euler can behave poorly for stiff choices of \( R \) and \( C \).  
- Improved Euler is stable for moderate step sizes.  
- Implicit Euler and Crank–Nicolson give smooth, stable curves.  
- Crank–Nicolson consistently provides the most accurate evolution.

Plots of capacitor voltage vs. time were produced for comparison.
---

## 3. HW 18 — Automatic Differentiation (AD)

In this part we implemented a minimal forward-mode automatic differentiation framework based on a templated class  
**`AutoDiff<N,T>`**, capable of tracking derivatives with respect to **N** variables.

### 3.1 Features Implemented

The AutoDiff class supports:
- A stored **value** `val`
- A stored **derivative vector** `deriv[0…N-1]`
- Initialization from:
  - constants
  - variables of the form `Variable<I>`
- Operator overloading for:
  - `+` addition  
  - `-` subtraction  
  - `*` multiplication (all 3 overloads: AD*AD, const*AD, AD*const)  
  - `/` division (full quotient rule support)
- Elementary functions:
  - `sin`, `cos`, `exp`, `log`, `pow`

This enables symbolic-like differentiation while evaluating normal C++ expressions.

### 3.2 Test Function

We tested the AD system on:

$$
f(x) = x^2 + 3\sin(x)
$$

Evaluated at \(x = 1\), the AD derivative was compared to a numeric finite-difference estimate.

#### **Result**
- AutoDiff derivative: **matches numerical derivative exactly**  
- Confirms correct implementation of forward-mode AD

### 3.3 Pendulum System 
A second AD test was implemented for the simple pendulum derivative:

$$
f(\theta) = -\frac{g}{L}\sin(\theta)
$$

AD successfully produced the correct symbolic derivative:

$$
f'(\theta) = -\frac{g}{L}\cos(\theta)
$$

and matched the numeric approximation.

---

## 4. HW 19 — Runge–Kutta Methods (RK2)

In this exercise we implemented **one explicit Runge–Kutta method** (RK2), as required by the assignment.  
RK4 was only used for internal testing and was removed before merging.

### 4.1 RK2 (Midpoint Method)
The RK2 method improves accuracy significantly compared to Euler-type methods and is defined by:

- Predictor step:  
  $$k_1 = f(y_n)$$

- Midpoint evaluation:  
  $$k_2 = f\left(y_n + \frac{1}{2}\tau k_1\right)$$

- Update:  
  $$y_{n+1} = y_n + \tau k_2$$

### 4.2 Testing

RK2 was tested on the mass–spring system for various step counts:

- 50 steps  
- 100 steps  
- 200 steps  
- 400 steps  

### 4.3 Observations

- **RK2** shows much smaller phase drift than Euler or improved Euler  
- The numerical trajectory becomes smoother as step size decreases  
- The energy error is significantly smaller than the Euler family  
- RK2 remains stable for much larger time steps than explicit Euler  

Overall, RK2 provides a clear accuracy improvement with minimal additional computation.

---

## 5. Overall Conclusions

- **Implicit Euler** is stable and suitable for stiff systems  
- **Crank–Nicolson** provides symmetric, high-quality solutions for oscillatory systems  
- **Automatic Differentiation** works correctly and matches finite-difference tests  
- **RK2** offers superior accuracy compared to Euler-type methods with only a small cost increase  
- The unified ODE framework allowed consistent testing across all numerical methods  



---
