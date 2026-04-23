import numpy as np
from params import Params

# -----------------------------
# grid generation
# -----------------------------
def grid_generation(p: Params):
    """
    Generate 1D grid and nozzle geometry.

    Returns:
        x    : grid coordinates
        A    : area distribution
        dAdx : area gradient
    """
    x = np.linspace(p.x0, p.x1, p.Imax + 1)

    dAdx = np.where(x <= 0.0, -0.0318, 0.0510)
    A = p.A_throat + dAdx * x

    return x, A, dAdx