import numpy as np

class NozzleGeom:
    """
    1D nozzle geometry and grid
    """

    def __init__(self, params, domain ):
        x = np.linspace(params.x0, params.x1, params.Imax + 1)

        self.x = x[domain.slice]

        # piecewise linear area gradient
        self.dAdx = np.where(self.x <= 0.0, -0.0318, 0.0510)

        # area distribution
        self.A = params.A_throat + self.dAdx * self.x