# Source-Vortex Panels Method & Momentum Integral based Methods to solve Airfoil

The main goal of this project is to conduct a quick airfoil analysis that can be used in preliminary design, first a potential flow problem is solved using source-vortex panels method to obtain lift coefficient and moment coefficient around quarter chord (c/4), then, Pohlhausen method is applied using the resulting external velocity from panels method to solve laminar boundary layer problem, head's method is then applied upon transition to turbulent flow.

## Dependencies
- Octave/Matlab.
- XFoil is used by default to generate the panels distribution, make sure it's path is defined in the terminal, to make sure it's callable from the terminal

#### Important Note:
If you're using windows, you might need to change the name of the xfoil program in the file "LOAD_AIRFOIL.m" from "xfoil" to "xfoil.exe"
 
## Project Files
### main.m
This is the main code that you should run, and it will run everything else automatically.
- You define the airfoil:
  - If it is NACA 4 digit, you can specify the number of panels nodes, and xfoil will generate the airfoil data into the folder ./MY_AIRFOILS/
  - Otherwise you'll have to put the .dat file in ./MY_AIRFOILS/ and provide
- Call the function SVPM in "SVPM.m" to solve panels method
- Call the function pohlhausen in pohlhausen.m to solve boundary layer problem for both the upper and lower surface
- Integrate the local friction coefficient to find total drag coefficient
- For the two problems, you can provide the number of panels nodes by setting the variable N (it is defined twice, once for panels method and once for pohlhausen) as an array so that you can check **Grid Convergence**

### LOAD_AIRFOIL.m
- This file writes some xfoil instructions into a file, passes it to xfoil, then xfoil generates the airfoil data and saves it as: "./MY_AIRFOILS/airfoil_name.dat"
- Matlab then reads the file, and you (optionally) chooses to close the trailing edge of the airfoil (by default, it doesn't)

### SVPM.m
Solves source vortex panels method, it solves a linear system of N equations, where N-1 is the number of panels on which a source is assumed to exist, and the last N-th equation comes from assuming a vortex with constant strength exists, also the last equation is the one used to impose the Kutta condition.

### boundary_layer.m
- Solves Laminar Boundary Layer problem using pohlhausen numerical method, it calculates the values of B.L. thickness, displacement thickness, momentum thickness, shear stress, and local friction coefficient which can be then integrated over distance to find the total drag coefficient
- Predicts the **Point of Transition**, and then switchs to turbulent boundary layer solving head's method using a 4th order runge-kutta scheme
- Predict **Point of Separation** (either laminar or turbulent) and then terminate the code

## RK4.m
- Runge Kutta Scheme

## f_dot.m
- The function that returns the derivatives of the state variables in Head's method, which are (theta and $Ue*theta*H$)

## G.m
- Function returns H1 from H

## HofH1.m
- The inversion of the H1=G(H) relation

## F.m
- A part of the differential equations that need to be solved in Head's method needs F(H1)

## Issues & Possible Future Work
- For boundary layer, some Results are not right, nevertheless, different from XFOIL.
- Results at high angle of attack are so wrong.
- Changing the velocity or chord length from unity causes problems, usually this needs to look carefully at normalization.
- Implement the xfoil algorithm of discretizing the airfoil, so that xfoil is no longer needed.
- Different panels distributions give very different results.

