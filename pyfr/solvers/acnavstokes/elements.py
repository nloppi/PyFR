# -*- coding: utf-8 -*-

from pyfr.solvers.aceuler.elements import BaseACFluidElements
from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.nputil import npeval


class ACNavierStokesElements(BaseACFluidElements,
                             BaseAdvectionDiffusionElements):
    def set_backend(self, backend, nscalupts, nonce):
        super().set_backend(backend, nscalupts, nonce)
        backend.pointwise.register('pyfr.solvers.acnavstokes.kernels.tflux')

        tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                       c=self.cfg.items_as('constants', float))

        sponge_type = self.cfg.get('sponge', 'sponge-type', 'none')

        if sponge_type == 'visc':
            tplargs['spng_vis'] = 'mu'
            spngmu_expr = self.cfg.get('sponge', 'sponge-mu')
            ploc_u = self.ploc_at_np('upts').swapaxes(0, 1)
            mu_u = npeval(spngmu_expr,
                {d: ploc_u[i] for i, d in enumerate('xyz'[:self.ndims])}
            )
            spngmu_upts = self._be.const_matrix(mu_u, tags={'align'})
            if 'flux' in self.antialias:
                ploc_q = self.ploc_at_np('qpts').swapaxes(0, 1)
                mu_q = npeval(spngmu_expr,
                    {d: ploc_q[i] for i, d in enumerate('xyz'[:self.ndims])}
                )
                spngmu_qpts = self._be.const_matrix(mu_q, tags={'align'})
            else:
                spngmu_qpts = None

        elif sponge_type == 'none':
            tplargs['spng_vis'] = 'none'
            spngmu_upts = spngmu_qpts = None

        else:
            raise ValueError('Invalid sponge-type option')

        if 'flux' in self.antialias:
            self.kernels['tdisf'] = lambda: backend.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, self.neles],
                u=self._scal_qpts, smats=self.smat_at('qpts'),
                f=self._vect_qpts, spngmu=spngmu_qpts
            )
        else:
            self.kernels['tdisf'] = lambda: backend.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, self.neles],
                u=self.scal_upts_inb, smats=self.smat_at('upts'),
                f=self._vect_upts, spngmu=spngmu_upts
            )
