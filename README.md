# MMA algorithm in FORTRAN

These subroutines are essentially Fortran adaptations/translations of the code initially developed by Professor Svanberg in Matlab (http://www.smoptit.se). The model parameters are specified in the initial lines of the MMA_Main.f90 and MMA_variables.f90 files.

Note: To execute the code, it is imperative to have the Lapack library installed, as the subroutines rely on the SGESV and DGESV functions. Additionally, a makefile is provided to streamline the compilation and execution process.

- Svanberg, K. (2002). A class of globally convergent optimization methods based on conservative convex separable approximations. SIAM J. Optim. 12 (2), p. 555-573
- Svanberg, K. (1987). The method of moving asymptotes—a new method for structural optimization. International journal for numerical methods in engineering, 24(2), 359-373.
