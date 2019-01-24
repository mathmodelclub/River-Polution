# River-Polution

## 1D Diffusion modeling

### Introduction
In this project we intend to model the spread of a chemical *u* in time and space along a one dimensional river. 
Imagine a one dimensional river, from *x = 0* to *x = L*. At a certain position *x = a* there is a chemical plant, which has a leakage. 
From there the pesticide will spread along the river, until it reaches the ocean.

                                  How is the substance *u(x,t)* going to evolve along space and time?

### Formulation
Mathematically this is a partial differential equation of mixed type, i.e. **parabolic**, maping the diffusion and **hyperbolic** maping the advection. 
The equation reads as: 

Where 
