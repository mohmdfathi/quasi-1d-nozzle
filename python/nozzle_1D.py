import numpy as np
from mpi4py import MPI
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
 
    comm = domain.comm
    rank = domain.rank
    size = domain.size

    for iter in range(pars.iter_max):

        # ---- CFL (local max + global reduction)
        wave_speed = np.abs(flow.u) + flow.sound_speed()
        dt_local = pars.CFL * pars.dx / np.max(wave_speed)
        dt = comm.allreduce(dt_local, op=MPI.MIN)
        dx = pars.dx
    
        # ---- conservative update (local)
        flow.update_conservatives(geom, Q)
        flow.update_fluxes_source(geom, F, H)
    
        # MacCormack predictor step (backward difference)
        dQdt_pred = (F[:, 1:-1] - F[:, 0:-2]) / dx + H[:, 1:-1]
        Qn[:, 1:-1] = Q[:, 1:-1] - dt * dQdt_pred
        add_dissipation(flow, Q, Qn)

        flow.update_primitives(geom, Qn) 
        flow.update_fluxes_source(geom, F, H)

        # MacCormack corrector step (forward difference)
        dQdt_corr = (F[:, 2:] - F[:, 1:-1]) / dx + H[:, 1:-1]
        dQdt_mean = 0.5 * (dQdt_pred + dQdt_corr)
        Qn[:, 1:-1] = Q[:, 1:-1] - dt * dQdt_mean
        add_dissipation(flow, Q, Qn)

        flow.update_primitives(geom, Qn) 

        if iter % 1000 == 0:

            if rank == size-1:
                comm.send(Qn[0, -2], dest=0, tag=11)

            if rank == 0:
                q_left_global = Qn[0, 1]
                q_right_global = comm.recv(source=size-1, tag=11)
                mdot_ref = pars.p0 * pars.A_throat * np.sqrt( ( pars.k / (pars.R * pars.T0)) * (2.0 / (pars.k + 1.0)) ** ((pars.k + 1.0) / (pars.k - 1.0)) )
                imbalance = (q_right_global - q_left_global) / mdot_ref
            
                print(f"[iter {iter}] mass imbalance = {imbalance:.4E}")

    filename = f"nozzle_maccormack_parallel_rank_{rank}.csv"
    write_results(flow, geom, filename)
