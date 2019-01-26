# River-Polution

## 1D Diffusion modeling

### Introduction
In this project we intend to model the spread of a chemical *u* in time and space along a one dimensional river. 
Imagine a one dimensional river, from *x = 0* to *x = L*. At a certain position *x = a* there is a chemical plant, which has a leakage. 
From there the pesticide will spread along the river, until it reaches the ocean.

                                  How is the substance u(x,t) going to evolve along space and time?

### Formulation
Mathematically this is a partial differential equation of mixed type, i.e. **parabolic**, maping the diffusion and **hyperbolic** maping the advection. 
The equation reads as: 
    <p align="center"> 
   <img src ="https://github.com/mathmodelclub/River-Polution/blob/Developer/Pics/Equation1.png" width="500" height="100" />    
   </p>


Where *v* is the velocity of the river, *D* is the diffusion coefficient, and *psi* is the source, i.e. the chemical diffusing into the river. 

### Modeling
We are using finite differences of second order accuracy to approximate the differetiation in space. This will lead to a tridiagonal square matr: 

   <p align="center"> 
   <img src ="https://github.com/mathmodelclub/River-Polution/blob/Developer/Pics/Matrix1.png" width="400" height="400" />    
   </p>

In the next step, we are using the Chrank Nicolson method (of order 2 in space and time) to mnumerically integrate. After doing so we find the following solution for the evolution of the pesticide in space and time.   



  Time animation of pesticide flow |  Contour depiciton of pesticide
  :-------------------------:|:-------------------------:
 <img src="https://im.ezgif.com/tmp/ezgif-1-5ad29534ab6f.gif" width="320" height="300" />  |  <img src="https://github.com/mathmodelclub/River-Polution/blob/Developer/Pics/Cont.PNG" width="320" height="300" />
 
 ### Environmental Impact
 We are analyzing the impact to the environment in the case the chemical plant leaks toxides along the river. In particular, we are interested in the number of dead fish along the river. Two different measures are taken into account:
 
 Total pesticide concentration | Maximal pesticide concentration
 
  
  
  
  
  
  
