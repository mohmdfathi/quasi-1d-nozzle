import numpy as np
from dataclasses import dataclass

@dataclass(frozen=True)
class Parameters:
    # -----------------------
    # solver settings
    # -----------------------
    Imax      : int = 1280
    iter_max  : int = 500_000

    s_shape   : tuple = (Imax+1,)
    v_shape   : tuple = (3,Imax+1)

    # -----------------------
    # physical constants
    # -----------------------
    R         : float = 286.9865
    k         : float = 1.4
    cv        : float = 718.0

    p0        : float = 25000.0
    T0        : float = 300.0
    pb        : float = 15000.0

    # -----------------------
    # numerical
    # -----------------------
    CFL       : float = 0.4
    Cx        : float = 5.0

    # -----------------------
    # geometry
    # -----------------------
    x0        : float =-0.25
    x1        : float = 1.03
    dx        : float = ( x1 - x0 ) / Imax
    A_throat  : float = 0.03150