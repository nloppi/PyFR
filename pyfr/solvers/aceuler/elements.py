# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionElements


class BaseFluidElements(object):
    formulations = ['std', 'dual']

    privarmap = {2: ['p', 'u', 'v'],
                 3: ['p', 'u', 'v', 'w']}

    convarmap = {2: ['p', 'u', 'v'],
                 3: ['p', 'u', 'v', 'w']}

    dualcoeffs = {2: ['u', 'v'],
                  3: ['u', 'v', 'w']}

    visvarmap = {
        2: {'velocity': ['u', 'v'],
            'pressure': ['p']},
        3: {'velocity': ['u', 'v', 'w'],
            'pressure': ['p']}
    }

    @staticmethod
    def pri_to_con(pris, cfg):
        return pris

    @staticmethod
    def con_to_pri(convs, cfg):
        return convs


class ACEulerElements(BaseFluidElements, BaseAdvectionElements):
    def set_backend(self, backend, nscalupts):
        super().set_backend(backend, nscalupts)

        # Register our flux kernel
        backend.pointwise.register('pyfr.solvers.aceuler.kernels.tflux')

        # Template parameters for the flux kernel
        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       c=self.cfg.items_as('constants', float))

        if 'flux' in self.antialias:
            self.kernels['tdisf'] = lambda: backend.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts
            )
        else:
            self.kernels['tdisf'] = lambda: backend.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts
            )
