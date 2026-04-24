import numpy as np
from params import Parameters 
from state  import FlowState
from grid   import NozzleGeom
from mpi    import Domain1D

from routines import write_results, add_dissipation

if __name__ == "__main__":

    pars = Parameters()
    domain = Domain1D(pars)
    flow = FlowState(pars, domain)
    geom = NozzleGeom(pars, domain)

    shape = domain.v_shape
    Q  = np.empty( shape )
    Qn = np.empty( shape )
    F  = np.empty( shape )
    H  = np.empty( shape )

    flow.initialize(mach=0.5)
 
    for iter in range(pars.iter_max):

        # ---- CFL time step
        wave_speed = np.abs(flow.u) + flow.sound_speed()
        dt = pars.CFL * pars.dx / np.max(wave_speed)
        dx = pars.dx
    
        # ---- build conservative forms
        flow.update_conservatives(geom, Q)
        flow.update_fluxes_source(geom, F, H)
    
        # MacCormack predictor step (backward difference)
        dQdt_pred = (F[:, 1:-1] - F[:, 0:-2]) / dx + H[:, 1:-1]
        Qn[:, 1:-1] = Q[:, 1:-1] - dt * dQdt_pred
        add_dissipation(flow, Qn)

        flow.update_primitives(geom, Qn) 
        flow.update_fluxes_source(geom, F, H)

        # MacCormack corrector step (forward difference)
        dQdt_corr = (F[:, 2:] - F[:, 1:-1]) / dx + H[:, 1:-1]
        dQdt_mean = 0.5 * (dQdt_pred + dQdt_corr)
        Qn[:, 1:-1] = Q[:, 1:-1] - dt * dQdt_mean
        add_dissipation(flow, Qn)

        flow.update_primitives(geom, Qn) 

        if iter % 1000 == 0:
            mdot_ref = pars.p0 * pars.A_throat * np.sqrt( ( pars.k / (pars.R * pars.T0)) * (2.0 / (pars.k + 1.0)) ** ((pars.k + 1.0) / (pars.k - 1.0)) )
            imbalance = (Qn[0, -2] - Qn[0, 0]) / mdot_ref
            print(f"[iter {iter}] mass imbalance = {imbalance:.4E}")

    write_results(flow, geom, "nozzle_maccormack_serial.csv")
