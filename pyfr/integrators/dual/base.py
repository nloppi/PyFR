# -*- coding: utf-8 -*-

from abc import abstractmethod

from pyfr.integrators.base import BaseIntegrator


class BaseDualIntegrator(BaseIntegrator):
    formulation = 'dual'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.dtaumin = 1.0e-12
        if self.cfg.get('solver-time-integrator', 'controller') == 'mg':
            self._dtaus = self._init_dtaus()
        else:
            self._dtau = self.cfg.getfloat('solver-time-integrator', 'pseudo-dt')

    def _init_dtaus(self):
        dtaus = self.cfg.getfloatlist('solver-time-integrator', 'pseudo-dts')
        return dtaus

    @property
    def _stepper_regidx(self):
        return self._regidx[0][:self._pseudo_stepper_nregs]

    @property
    def _source_regidx(self):
        return self._regidx[0][self._pseudo_stepper_nregs:]

    @property
    def _MG_regidx(self):
        return self._regidx[0][-3:]

    @abstractmethod
    def _dual_time_source(self):
        pass

    @abstractmethod
    def finalise_step(self, currsoln):
        pass
