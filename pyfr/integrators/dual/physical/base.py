# -*- coding: utf-8 -*-

from pyfr.integrators.base import BaseIntegrator
from pyfr.integrators.dual.pseudo import get_pseudo_integrator
from pyfr.util import proxylist


class BaseDualIntegrator(BaseIntegrator):
    formulation = 'dual'

    def __init__(self, backend, systemcls, rallocs, mesh, initsoln, cfg):
        # Composite the base classes together to form a new type
        self.pseudointegrator = get_pseudo_integrator(
            backend, systemcls, rallocs, mesh,
            initsoln, cfg, self._stepper_coeffs
        )

        super().__init__(backend, rallocs, mesh, initsoln, cfg)

        # Event handlers for advance_to
        self.completed_step_handlers = proxylist(
            self._get_plugins()
        )
        
        # Delete the memory-intensive elements map from the system
        del self.system.ele_map

    @property
    def system(self):
        return self.pseudointegrator.system

    @property
    def _stepper_coeffs(self):
        pass

    @property
    def pseudostepinfo(self):
        # Stats on the most recent step
        return self.pseudointegrator.pseudostepinfo

    @property
    def soln(self):
        # If we do not have the solution cached then fetch it
        if not self._curr_soln:
            self._curr_soln = self.pseudointegrator.system.ele_scal_upts(
                self.pseudointegrator._idxcurr
            )

        return self._curr_soln
