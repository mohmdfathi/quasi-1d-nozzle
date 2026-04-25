# quasi-1d-nozzle

A nozzle is a shaped duct that converts pressure energy into kinetic energy by changing the cross-sectional area of the flow passage. Nozzles are essential in many engineering systems because they accelerate fluids, control mass flow rate, generate thrust, and enable supersonic motion. Convergent-divergent nozzles are especially important in applications such as rocket engines, gas turbines, propulsion devices, and wind tunnels.

The flow inside a convergent-divergent nozzle strongly depends on the imposed back pressure. Under specific outlet pressure conditions, the flow can accelerate to supersonic speeds and a normal shock may form inside the divergent section, causing an abrupt transition back to subsonic flow downstream. Capturing such discontinuities accurately is one of the main challenges in the numerical solution of hyperbolic conservation laws, even when the initial conditions are smooth.

This project numerically solves the quasi-1D compressible Euler equations for flow through a convergent-divergent nozzle using the explicit MacCormack scheme. The solver is implemented in Python, Fortran, and C++, with  simple plans for CPU and GPU parallelization.

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