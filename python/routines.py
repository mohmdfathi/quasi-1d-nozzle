"""
==============================================================
MODULE: routines
==============================================================

This module implements the core numerical routines for a
quasi-1D compressible Euler solver.

--------------------------------------------------------------
CONTENTS
--------------------------------------------------------------

1. primitives_to_conservatives
   Convert primitive variables (rho, u, p) → conservative variables Q

2. conservatives_to_primitives
   Convert conservative variables Q → primitive variables (rho, u, p)

3. apply_boundaries
   Enforce inlet/outlet boundary conditions:
   Subsonic inlet and outlet pressure

4. fluxes_and_sources
   Compute: Convective flux vector F & Geometric source term H

5. add_dissipation
   Jameson-type artificial dissipation for shock capturing
   and numerical stabilization

6. write_results
   Post-processing routine:
   - nondimensionalization
   - diagnostic quantities (Mach, mass flow, pressure ratio)
   - export to ASCII (Fortran-compatible format)

==============================================================
"""

import numpy as np

# ==============================================================
def primitives_to_conservatives(flow, geom, Q):
    """
    Convert primitive variables → conservative variables
    """

    k = flow.params.k

    rho = flow.rho
    u   = flow.u
    p   = flow.p
    A   = geom.A

    E = p / ((k - 1.0) * rho) + 0.5 * u**2

    Q[0] = rho * A
    Q[1] = rho * u * A
    Q[2] = rho * E * A


# ==============================================================
def apply_boundaries(flow):
    """
    Apply inlet/outlet BCs (in-place)
    """

    rho = flow.rho
    u   = flow.u
    p   = flow.p

    # [ ghost | interior ... interior | ghost ]
    #   0        1 ... N-2              N-1
    # neighbors
    left  = flow.domain.left
    right = flow.domain.right

    comm  = flow.domain.comm


    # pack
    flow.send_left[:]  = (rho[1],  u[1],  p[1])
    flow.send_right[:] = (rho[-2], u[-2], p[-2])

    # exchange
    comm.Sendrecv(flow.send_left,  left,  recvbuf=flow.recv_left,  source=left)
    comm.Sendrecv(flow.send_right, right, recvbuf=flow.recv_right, source=right)

    # unpack
    rho[0],  u[0],  p[0]  = flow.recv_left
    rho[-1], u[-1], p[-1] = flow.recv_right


    par = flow.params
    p0, T0, pb, R, k = par.p0, par.T0, par.pb, par.R, par.k

    # inlet
    if( flow.domain.rank == 0 ):
        p[0] = p[1]

        M2 = ((p0 / p[0])**((k - 1.0) / k) - 1.0) * (2.0 / (k - 1.0))
        T  = T0 / (1.0 + 0.5 * (k - 1.0) * M2)

        rho[0] = p[0] / (R * T)
        u[0]   = np.sqrt(M2 * k * R * T)

    # outlet
    if( flow.domain.rank == flow.domain.size - 1 ):
        p[-1]   = pb
        rho[-1] = 2.0 * rho[-2] - rho[-3]
        u[-1]   = 2.0 * u[-2]   - u[-3]


# ==============================================================
def conservatives_to_primitives(flow, geom, Q):
    """
    Update primitive variables only on interior nodes (1:-1)
    Boundaries are handled separately via BCs
    """

    k = flow.params.k
    A = geom.A

    # interior slice
    i = slice(1, -1)

    flow.rho[i] = Q[0, i] / A[i]
    flow.u[i]   = Q[1, i] / Q[0, i]

    E = Q[2, i] / Q[0, i]

    flow.p[i] = (k - 1.0) * flow.rho[i] * (E - 0.5 * flow.u[i]**2)

    apply_boundaries(flow)


# ==============================================================
def fluxes_and_sources(flow, geom, F, H):
    """
    Compute fluxes and source terms
    """

    k = flow.params.k

    rho = flow.rho
    u   = flow.u
    p   = flow.p

    A    = geom.A
    dAdx = geom.dAdx

    E = p / ((k - 1.0) * rho) + 0.5 * u**2

    F[0] = rho * u * A
    F[1] = (rho * u**2 + p) * A
    F[2] = u * (rho * E + p) * A

    H[0] = 0.0
    H[1] = -p * dAdx
    H[2] = 0.0


# ==============================================================
def add_dissipation(flow, Qn):
    """
    Jameson artificial dissipation (restricted interior stencil)
    """

    C_visc = flow.params.Cx
    p = flow.p

    i = slice(2, -2)

    nu = np.abs( p[2:] - 2*p[1:-1] + p[:-2] ) / \
               ( p[2:] + 2*p[1:-1] + p[:-2] )

    D2 = Qn[:, 2:] - 2*Qn[:, 1:-1] + Qn[:, :-2]

    Qn[:, i] += C_visc * nu[1:-1] * D2[:,1:-1]
    

# ==============================================================
def write_results(flow, geom, filename):
    
    par = flow.params
    p0, T0, R, k     = par.p0, par.T0, par.R, par.k
    A_throat, x0, x1 = par.A_throat, par.x0, par.x1

    rho, u, p = flow.rho, flow.u, flow.p

    x, A = geom.x, geom.A

    # =========================
    # reference values
    # =========================
    rho_0 = p0 / (R * T0)
    a_0 = np.sqrt(k * R * T0)
    mdot_ref = rho_0 * a_0 * A_throat

    # =========================
    # vectorized calculations
    # =========================
    x_over_L = x / (x1 - x0)
    p_over_p0 = p / p0
    Mach = u / np.sqrt(k * p / rho)
    mdot_nd = (rho * u * A) / mdot_ref

    # =========================
    # stack like Fortran columns
    # =========================
    data = np.column_stack( (x_over_L, p_over_p0, Mach, mdot_nd, A) )

    # =========================
    # write file (Fortran-style simple write)
    # =========================
    with open(filename, "w") as f:
        f.write("x/L,p/p0,Mach,mdot/mdot_ref,A\n")
        np.savetxt(f, data, fmt=["%16.8E"] * data.shape[1])