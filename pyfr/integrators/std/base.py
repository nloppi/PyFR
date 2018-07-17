# -*- coding: utf-8 -*-

from pyfr.integrators.base import BaseIntegrator
from pyfr.integrators.base import BaseCommon
from pyfr.util import proxylist


class BaseStdIntegrator(BaseCommon, BaseIntegrator):
    formulation = 'std'

    def __init__(self, backend, systemcls, rallocs, mesh, initsoln, cfg):
        super().__init__(backend, rallocs, mesh, initsoln, cfg)
        # Determine the amount of temp storage required by this method
        self.nregs = self._stepper_nregs

        # Construct the relevant mesh partition
        self.system = systemcls(backend, rallocs, mesh, initsoln,
                                nregs=self.nregs, cfg=cfg)

        # Storage for register banks and current index
        self._init_reg_banks(backend, self.system)

        # Global degree of freedom count
        self._gndofs = self._get_gndofs(self.system)

        # Add kernel cache
        self._axnpby_kerns = {}

        # Event handlers for advance_to
        self.completed_step_handlers = proxylist(self._get_plugins())

        # Sanity checks
        if self._controller_needs_errest and not self._stepper_has_errest:
            raise TypeError('Incompatible stepper/controller combination')

        # Ensure the system is compatible with our formulation
        if self.formulation not in systemcls.elementscls.formulations:
            raise RuntimeError(
                'System {0} does not support time stepping formulation {1}'
                .format(systemcls.name, self.formulation)
            )

    @property
    def soln(self):
        # If we do not have the solution cached then fetch it
        if not self._curr_soln:
            self._curr_soln = self.system.ele_scal_upts(self._idxcurr)

        return self._curr_soln

    @property
    def _controller_needs_errest(self):
        pass

    @property
    def _stepper_has_errest(self):
        pass
