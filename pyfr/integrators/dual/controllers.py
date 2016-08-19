# -*- coding: utf-8 -*-

import math

from pyfr.integrators.dual.base import BaseDualIntegrator
from pyfr.mpiutil import get_comm_rank_root, get_mpi
from pyfr.util import memoize


class BaseDualController(BaseDualIntegrator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Solution filtering frequency
        self._fnsteps = self.cfg.getint('soln-filter', 'nsteps', '0')

        # Stats on the most recent step
        self.stepinfo = []

    def _accept_step(self, dt, idxcurr):
        self.tcurr += dt
        self.nacptsteps += 1
        self.nacptchain += 1

        self._idxcurr = idxcurr

        # Filter
        if self._fnsteps and self.nacptsteps % self._fnsteps == 0:
            self.system.filt(idxcurr)

        # Invalidate the solution cache
        self._curr_soln = None

        # Fire off any event handlers
        self.completed_step_handlers(self)

    @property
    def nsteps(self):
        return self.nacptsteps + self.nrjctsteps


class DualNoneController(BaseDualController):
    controller_name = 'none'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        sect = 'solver-time-integrator'

        self._innersteps = 0
        self._idxprev = 0
        self._maxniters = self.cfg.getint(sect, 'niters')
        self._atol = self.cfg.getfloat(sect, 'atol', default=1e-6)
        self._rtol = self.cfg.getfloat(sect, 'rtol', default=1e-6)
        self._chkevery = self.cfg.getfloat(sect, 'chk', default=5)

    def advance_to(self, t):
        if t < self.tcurr:
            raise ValueError('Advance time is in the past')

        while self.tcurr < t:
            for i in range(self._maxniters):
                dt = max(min(t - self.tcurr, self._dt), self.dtmin)
                dtau = max(min(t - self.tcurr, self._dtau), self.dtaumin)

                # Take the step
                self._idxcurr, self._idxprev = self.step(self.tcurr, dt, dtau)
                self._innersteps += 1
                if self._innersteps % self._chkevery == 0:
                    # Subtract the current solution from the previous solution
                    self._add(1.0, self._idxprev, -1.0, self._idxcurr)
                    # Compute the normalised residual
                    resid = self._resid(self._idxprev, self._idxcurr)
                    if resid < 1.0:
                        break

            # Update the dual-time stepping banks (n+1 => n, n => n-1)
            self.finalise_step(self._idxcurr)

            # We are not adaptive, so accept every step
            self._accept_step(dt, self._idxcurr)

    @memoize
    def _get_errest_kerns(self):
        return self._get_kernels('errest', nargs=3, norm='uniform')

    def _resid(self, x, y):
        comm, rank, root = get_comm_rank_root()

        # Get an errest kern to compute the maximum residual
        errest = self._get_errest_kerns()

        # Prepare and run the kernel
        self._prepare_reg_banks(x, y, y)
        self._queue % errest(self._atol, self._rtol)

        # Reduce locally (element types) and globally (MPI ranks)
        rl = max(errest.retval)
        rg = comm.allreduce(rl, op=get_mpi('max'))

        # Normalise
        return math.sqrt(rg)
