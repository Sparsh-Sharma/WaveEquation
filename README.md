# Solving the Wave Equation
This package is designed to solve the wave equation

![](https://latex.codecogs.com/svg.latex?\frac{\partial^2%20A}{\partial%20t^2}=c^2\frac{\partial^2%20A}{\partial%20x^2})

## Assignment
Create a Julia package that solves the scalar wave equation in 1+1 dimensions. You will need to:

- choose a domain (e.g. the unit interval [0, 1])
- choose a suitable discretization method (e.g. piecewise linear continuous functions) for this domain
- choose a formulation of the wave equation, i.e. a concrete set of evolved variable ("state vector") and evolution equations for these
- define a suitable (conserved) energy
- choose suitable initial and boundary conditions
- implement this in Julia
- compare solutions at different resolutions and demonstrate convergence
- create a figure(s) so that others can understand your work
- present all the above in a git repository

## Method
The wave equation is solved on a domain of [0,1] using simple position discretization (piecewise linear continuous functions). An initial configuration and a first time derivative for the wavefunction *A* are specified

The wave equation can be represented
d/dt (A, dA/dt) = (dA/dt, d2A/dx2)

We use 2nd order Runge-Kutta integration to step forward this equation in time
which is a formalism agnostic of the spatial discretization

d/dt (state) = f(state)
