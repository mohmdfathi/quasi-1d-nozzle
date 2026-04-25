# quasi-1d-nozzle
A simple simulation of quasi-1D compressible nozzle flow implemented in Python, Fortran, and C++, with plans for CPU and GPU parallelization.


## Governing Equations

We solve the quasi-1D unsteady compressible flow in ducts of variable cross-section. The governing equations are:

\[
\partial_t Q + \partial_x F + H = 0
\]

where the solution vector \( Q \) and flux vector \( F \) can be computed using the primitive variables \( \rho, u, p\), which are density, velocity, static pressure of the flow:

\[ Q = \begin{pmatrix} \rho A \\ \rho u A \\ \rho E A \end{pmatrix} ,\quad
   F = \begin{pmatrix} \rho u A \\(\rho u^2 + p) A \\u \left( \rho E + p \right) A \end{pmatrix} \]

Here \(A\) denotes the cross section area of the duct, and \(\gamma\) is the specific heat ratio of the gas. The source term accounts for the effect of the cross section variation on the inviscid flow:

\[ H = \begin{pmatrix} 0 \\ -p \frac{dA}{dx}\\ 0 \end{pmatrix} \]

This set of equations represents conservation of mass, momentum, and energy. To close the system, we use the equation of state:

\[ E = \frac{u^2}{2} + \frac{p}{\rho(\gamma - 1)} \]

where \( E \) is the total specific internal enegy of the gas including the kinetic energy.

### Boundary Conditions

At the inlet, total (stagnation) conditions are prescribed in terms of total pressure \(P_0\) and total temperature \(T_0\). These values are used to reconstruct the primitive flow variables \((\rho, u, p)\) assuming isentropic relations, with appropriate treatment of characteristic waves for subsonic inflow where only the incoming characteristic is imposed while the outgoing information is taken from the interior solution. At the outlet, a fixed static back pressure \(p_{b}\) is imposed. For subsonic outflow, only one characteristic is specified through the back pressure while the remaining variables are extrapolated from the interior, whereas for supersonic outflow all variables are fully extrapolated from the interior domain. 