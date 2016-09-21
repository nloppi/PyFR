# -*- coding: utf-8 -*-

from pyfr.integrators.dual.base import BaseDualIntegrator


class BaseDualStepper(BaseDualIntegrator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        elementscls = self.system[0].elementscls
        self._subdims = [elementscls.convarmap[self.system[0].ndims].index(v)
                         for v in elementscls.dualcoeffs[self.system[0].ndims]]

    @property
    def _stepper_nregs(self):
        return self._pseudo_stepper_nregs + len(self._dual_time_source) - 1


class BDF2DualStepper(BaseDualStepper):
    stepper_name = 'bdf2'

    @property
    def _stepper_order(self):
        return 2

    @property
    def _dual_time_source(self):
        return [-1.5, 2.0, -0.5]


class BDF3DualStepper(BaseDualStepper):
    stepper_name = 'bdf3'

    @property
    def _stepper_order(self):
        return 3

    @property
    def _dual_time_source(self):
        return [-11.0/6.0, 3.0, -1.5, 1.0/3.0]


class BackwardEulerDualStepper(BaseDualStepper):
    stepper_name = 'backward-euler'

    @property
    def _stepper_order(self):
        return 1

    @property
    def _dual_time_source(self):
        return [-1.0, 1.0]
