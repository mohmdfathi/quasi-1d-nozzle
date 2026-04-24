import numpy as np
from routines import *

class FlowState:
    """
    Container for primitive flow variables + parameters
    """

    def __init__(self, params, domain):
        self.params = params   # attach params
        self.domain = domain   # attach domain

        shape = self.domain.s_shape

        self.rho = np.empty(shape)
        self.u   = np.empty(shape)
        self.p   = np.empty(shape)

        self.send_left  = np.empty(3)
        self.recv_left  = np.empty(3)
        self.send_right = np.empty(3)
        self.recv_right = np.empty(3)

    def initialize(self, mach=0.5):
        """
        Initialize uniform isentropic flow
        """

        k = self.params.k
        R = self.params.R
        p0 = self.params.p0
        T0 = self.params.T0

        fac = 1.0 + 0.5 * (k - 1.0) * mach**2

        self.p[:] = p0 / fac**(k / (k - 1.0))
        T = T0 / fac

        self.rho[:] = self.p / (R * T)
        self.u[:] = mach * np.sqrt(k * R * T)

    def sound_speed(self):
        k = self.params.k
        return np.sqrt(k * self.p / self.rho)
    
    def update_conservatives(self, geom, Q):
        return primitives_to_conservatives(self, geom, Q)
    
    def update_fluxes_source(self, geom, F, H):
        return fluxes_and_sources(self, geom, F, H)
    
    def update_primitives(self, geom, Q):
        return conservatives_to_primitives(self, geom, Q)