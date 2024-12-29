# Source-Vortex Panels Method & Pohlhausen Method to solve Airfoil

The main goal of this project is to conduct a quick airfoil analysis that can be used in preliminary design, first a potential flow problem is solved using source-vortex panels method to obtain lift coefficient and moment coefficient around quarter chord (c/4), then, Pohlhausen method is applied using the resulting external velocity from panels method to solve boundary layer problem and determine drag.
Up until now, the code only calculates laminar drag, it also predicts the transition point and separation point, but **no** model for turbulence is made.

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

### velocity.m
This code is needed before solving the boundary layer problem, because B.L. problem needs exernal velocity to be known, the code interpolates the velocity values from panels method, then applies a **4th order central difference** to find the 1st and 2nd derivatives of velocity, which are needed by the pohlhausen method
- CD-4 is used to reduce the high instabilities that occur when the grid points increase, specially at high angles of attack

### pohlhausen.m
- Solves Boundary Layer problem using pohlhausen numerical method, it calculates the values of B.L. thickness, displacement thickness, momentum thickness, shear stress, and local friction coefficient which can be then integrated over distance to find the total drag coefficient
- Predicts the **Point of Transition**, but still continue to solve laminar flow
- Predict **Point of Separation** and then terminate the code

## Possible Future Work (ordered by priority)
- Separation seems to usually happen early in the program
- Implement the xfoil algorithm of making panels distribution, so that xfoil is no longer needed
- Different panels distributions give very different results
- Include code for solving turbulent boundary layer flow (Head's Method)
