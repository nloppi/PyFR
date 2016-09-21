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

        self._idxcurr[0] = idxcurr

        # Filter
        if self._fnsteps and self.nacptsteps % self._fnsteps == 0:
            self.system[0].filt(idxcurr)

        # Invalidate the solution cache
        self._curr_soln = None

        # Fire off any event handlers
        self.completed_step_handlers(self)


class DualNoneController(BaseDualController):
    controller_name = 'none'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        sect = 'solver-time-integrator'

        self._maxniters = self.cfg.getint(sect, 'pseudo-niters-max')
        self._minniters = self.cfg.getint(sect, 'pseudo-niters-min')

        if self._maxniters < self._minniters:
            raise ValueError('The maximum number of pseudo-iterations must '
                             'be greater than or equal to the minimum')

        self._pseudo_aresid = self.cfg.getfloat(sect, 'pseudo-aresid')
        self._pseudo_rresid = self.cfg.getfloat(sect, 'pseudo-rresid')

    def advance_to(self, t):
        if t < self.tcurr:
            raise ValueError('Advance time is in the past')

        while self.tcurr < t:
            for i in range(self._maxniters):
                dt = max(min(t - self.tcurr, self._dt), self.dtmin)
                dtau = max(min(t - self.tcurr, self._dtau), self.dtaumin)

                # Take the step
                self._idxcurr[0], self._idxprev = self.step(self.tcurr, dt, dtau)

                # Activate convergence monitoring after pseudo-niters-min
                if i >= self._minniters - 1:
                    # Subtract the current solution from the previous solution
                    self._add(-1.0/dtau, self._idxprev, 1.0/dtau, self._idxcurr[0])

                    # Compute the normalised residual and check for convergence
                    if self._resid(self._idxprev, self._idxcurr[0]) < 1.0:
                        break

            # Update the dual-time stepping banks (n+1 => n, n => n-1)
            self.finalise_step(self._idxcurr[0])

            # We are not adaptive, so accept every step
            self._accept_step(dt, self._idxcurr[0])

    @memoize
    def _get_errest_kerns(self):
        return self._get_kernels('errest', nargs=3, norm='uniform')

    def _resid(self, x, y):
        comm, rank, root = get_comm_rank_root()

        # Get an errest kern to compute the square of the maximum residual
        errest = self._get_errest_kerns()

        # Prepare and run the kernel
        self._prepare_reg_banks(x, y, y)
        self._queue % errest(self._pseudo_aresid, self._pseudo_rresid)

        # Reduce locally (element types) and globally (MPI ranks)
        rl = max(errest.retval)
        rg = comm.allreduce(rl, op=get_mpi('max'))

        # Normalise
        return math.sqrt(rg)


class DualMGController(BaseDualController):
    controller_name = 'mg'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._niters = self.cfg.getint('solver-time-integrator', 'mg-cycles')

    def advance_to(self, t):
        if t < self.tcurr:
            raise ValueError('Advance time is in the past')

        while self.tcurr < t:
            dt = max(min(t - self.tcurr, self._dt), self.dtmin)

            # Number of V-cycles
            for i in range(self._niters):

                # V-cycle down
                for j in range(self.levels - 1):
                    for k in range(self.leveliters[j]):
                        dtau = max(min(t - self.tcurr, self._dtaus[j]), self.dtaumin)
                        self._idxcurr[j], self._idxprev = self.step(self.tcurr, dt, dtau, level=j)
                    self.restrict(j, j+1)

                # V-cycle up
                for j in range(self.levels - 1):
                    for k in range(self.leveliters[self.levels - 1 - j]):
                        dtau = max(min(t - self.tcurr, self._dtaus[self.levels - 1 - j]), self.dtaumin)
                        self._idxcurr[self.levels - 1 - j], self._idxprev = self.step(self.tcurr, dt, dtau, level=(self.levels - 1 - j))
                    self.prolongate(self.levels - 1 - j, self.levels - 1 - (j + 1))

                # Post-smoothing at the highest level
                for k in range(self.leveliters[0]):
                    dtau = max(min(t - self.tcurr, self._dtaus[0]), self.dtaumin)
                    self._idxcurr[0], self._idxprev = self.step(self.tcurr, dt, dtau, level=0)

            # Update the dual-time stepping banks (n+1 => n, n => n-1)
            self.finalise_step(self._idxcurr[0])

            # We are not adaptive, so accept every step
            self._accept_step(dt, self._idxcurr[0])