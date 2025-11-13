# Welcome to ASC-ODE's documentation!


ASC-ODE is is a C++ library for solving ordinary differential equations.
The equation is defined by the right hand side function. ASC-ODE provides various time-steppers
which may be used as follos ...


```cpp
double tend = 4*M_PI;
int steps = 100;
double tau = tend/steps;

Vector<> y = { 1, 0 };  // initializer list
shared_ptr<NonlienarFunction> rhs = std::make_shared<MassSpring>(mass, stiffness);
  
ExplicitEuler stepper(rhs);

std::cout << 0.0 << "  " << y(0) << " " << y(1) << std::endl;
for (int i = 0; i < steps; i++)
  {
     stepper.DoStep(tau, y);
     std::cout << (i+1) * tau << "  " << y(0) << " " << y(1) << std::endl;
  }
```    

The result of the simulation in phase space is shown here

```{image} pictures/massspring_phase.png
:width: 40%
:align: center
```




## Installation

install it via git-clone:

    git clone https://github.com/my-github-clone/my-ode-solver.git


To configure and build some tests do

    cd my-ode-solver
    mkdir build
    cd build
    cmake ..
    make
    

## Available time-stepping methods are
...




   
