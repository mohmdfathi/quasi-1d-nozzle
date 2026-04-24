import numpy as np
from mpi4py import MPI

class Domain1D:
    def __init__(self, pars):

        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        # ---------- global info ----------
        self.N_global = pars.Imax + 1

        # ---------- decomposition ----------
        i_start, i_end = self._decompose( self.N_global )

        # ---------- neighbors ----------
        self.left  = self.rank - 1 if self.rank > 0 else MPI.PROC_NULL
        self.right = self.rank + 1 if self.rank < self.size - 1 else MPI.PROC_NULL

        # ---------- ghost settings ----------
        ng = 1
        self.ng_left  = 0 if self.rank == 0 else ng
        self.ng_right = 0 if self.rank == self.size - 1 else ng

        self.N_tot = i_end + self.ng_right - (i_start - self.ng_left)
        self.slice = slice(i_start-self.ng_left, i_end+self.ng_right)

        # ---------- shapes ----------
        nv = 3
        self.s_shape = (self.N_tot,)
        self.v_shape = (nv, self.N_tot)

    # -----------------------------
    def _decompose(self, N):

        nproc = self.size # number of processes
        rank  = self.rank # rank of the current process

        indices = np.array_split(np.arange(N), nproc)[rank]

        return indices[0], indices[-1] + 1