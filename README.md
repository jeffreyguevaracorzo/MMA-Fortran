# MMA Algorithm in FORTRAN

These subroutines are a Fortran adaptation of the code originally developed by Professor Svanberg in Matlab (http://www.smoptit.se). The model parameters are defined in the initial lines of the MMA_Main.f90 and MMA_variables.f90 files.

**Note:** To run the code, the LAPACK library must be installed, as the subroutines depend on the SGESV and DGESV functions. A makefile is included to facilitate the compilation and execution process. Feedback on code performance improvements or bug reports is highly appreciated. This code is intended for academic purposes only.

References:
- Svanberg, K. (2002). *A class of globally convergent optimization methods based on conservative convex separable approximations*. SIAM Journal on Optimization, 12(2), 555-573.
- Svanberg, K. (1987). *The method of moving asymptotes—a new method for structural optimization*. International Journal for Numerical Methods in Engineering, 24(2), 359-373.
