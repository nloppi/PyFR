# -*- coding: utf-8 -*-

from collections import defaultdict
import itertools as it
import re

from pyfr.inifile import Inifile
from pyfr.integrators.dual.pseudo.base import BaseDualPseudoIntegrator
from pyfr.util import memoize, proxylist


class DualMultiPIntegrator(BaseDualPseudoIntegrator):
    def __init__(self, backend, systemcls, rallocs, mesh, initsoln, cfg,
                 tcoeffs, multip_types):
        self.backend = backend

        sect = 'solver-dual-time-integrator-multip'
        # Get the solver order
        self._order = cfg.getint('solver', 'order')

        # Get the multigrid cycle
        self.cycle, self.csteps = zip(*cfg.getliteral(sect, 'cycle'))
        self.levels = sorted(set(self.cycle), reverse=True)
        self.level = self._order

        if max(self.cycle) > self._order:
            raise ValueError('The multigrid level orders cannot exceed '
                             'the solution order')

        if any(abs(i - j) > 1 for i, j in zip(self.cycle, self.cycle[1:])):
            raise ValueError('The orders of consecutive multigrid levels can '
                             'only change by one')

        if self.cycle[0] != self._order or self.cycle[-1] != self._order:
            raise ValueError('The multigrid cycle needs to start end with the '
                             'highest (solution) order ')

        # Multigrid pseudo-time steps
        dtau = cfg.getfloat('solver-time-integrator', 'pseudo-dt')
        dtauf = cfg.getfloat(sect, 'pseudo-dt-fact', 1.0)
        self.dtaus = {l: dtau*dtauf**(self._order - l) for l in self.levels}

        self._maxniters = cfg.getint('solver-time-integrator',
                                     'pseudo-niters-max', 0)
        self._minniters = cfg.getint('solver-time-integrator',
                                     'pseudo-niters-min', 0)

        # Generate suitable config files for lower multigrid levels
        mgcfgs = {l: Inifile(cfg.tostr()) for l in self.levels[1:]}
        for l, mgcfg in mgcfgs.items():
            mgcfg.set('solver', 'order', l)
            mgcfg.set('solver-time-integrator', 'pseudo-dt', self.dtaus[l])

            for sec in cfg.sections():
                m = re.match(r'solver-(.*)-mg-p{0}'.format(l), sec)
                if m:
                    mgcfg.rename_section(m.group(0), 'solver-' + m.group(1))

        # Insert the original config file to the multigrid config dictionary
        mgcfgs[self._order] = cfg

        # Get pseudo-integrators at all multigrid levels
        self.integrators = {}
        for i, psint in multip_types.items():
            # A class that bypasses pseudo-controller methods within a cycle
            class lpsint(psint):
                aux_nregs = 2 if i != self._order else 0

                @property
                def _aux_regidx(self):
                    return (self._regidx[-2:]
                            if self.aux_nregs is not 0 else None)

                def conv_mon(self, *args, **kwargs):
                    pass

                def finalise_pseudo_advance(self, *args, **kwargs):
                    pass

                def _rhs_with_dts(self, t, uin, fout, c=1):
                    # Compute -∇·f
                    self.system.rhs(t, uin, fout)

                    # Coefficients for the physical stepper
                    svals = [c*sc for sc in self._stepper_coeffs]

                    # Physical stepper source addition -∇·f - dQ/dt
                    axnpby = self._get_axnpby_kerns(len(svals) + 1,
                                                    subdims=self._subdims)
                    self._prepare_reg_banks(fout, self._idxcurr,
                                            *self._stepper_regidx)
                    self._queue % axnpby(1, *svals)

                    # Multigrid r addition
                    if self._aux_regidx:
                        axnpby = self._get_axnpby_kerns(2)
                        self._prepare_reg_banks(fout, self._aux_regidx[0])
                        self._queue % axnpby(1, -1)

            self.integrators[i] = lpsint(backend, systemcls, rallocs, mesh,
                                         initsoln, mgcfgs[i], tcoeffs)

        # Get a convergence monitoring method from the highest level controller
        self.mg_conv_mon = multip_types[self._order].conv_mon

        # Get the highest p system from plugins
        self.system = self.integrators[self._order].system

        # Initialise the restriction and prolongation matrices
        self._init_proj_mats()

        # Delete remaining elements maps from multigrid systems
        for l in self.levels[1:]:
            del self.integrators[l].system.ele_map

    def finalise_mg_advance(self, currsoln):
        psnregs = self.integrator._pseudo_stepper_nregs
        snregs = self.integrator._stepper_nregs

        # Rotate the stepper registers to the right by one
        self.integrator._regidx[psnregs:psnregs + snregs] = (
            self.integrator._stepper_regidx[-1:] +
            self.integrator._stepper_regidx[:-1]
        )

        # Copy the current soln into the first source register
        self.integrator._add(0, self.integrator._regidx[psnregs], 1, currsoln)

    @property
    def _idxcurr(self):
        return self.integrator._idxcurr

    @_idxcurr.setter
    def _idxcurr(self, y):
        self.integrator._idxcurr = y

    @property
    def pseudostepinfo(self):
        return self.integrator.pseudostepinfo

    @pseudostepinfo.setter
    def pseudostepinfo(self, y):
        self.integrator.pseudostepinfo = y

    @property
    def integrator(self):
        return self.integrators[self.level]

    def _init_proj_mats(self):
        self.projmats = defaultdict(proxylist)
        cmat = lambda m: self.backend.const_matrix(m, tags={'align'})

        for l in self.levels[1:]:
            for etype in self.integrator.system.ele_types:
                b1 = self.integrators[l].system.ele_map[etype].basis.ubasis
                b2 = self.integrators[l + 1].system.ele_map[etype].basis.ubasis

                self.projmats[l, l + 1].append(cmat(b1.proj_to(b2)))
                self.projmats[l + 1, l].append(cmat(b2.proj_to(b1)))

    @memoize
    def mgproject(self, l1, l2):
        inbanks = self.integrators[l1].system.eles_scal_upts_inb
        outbanks = self.integrators[l2].system.eles_scal_upts_inb

        return proxylist(
            self.backend.kernel('mul', pm, inb, out=outb)
            for pm, inb, outb in zip(self.projmats[l1, l2], inbanks, outbanks)
        )

    def restrict(self, l1, l2):
        l1idxcurr = self.integrators[l1]._idxcurr
        l2idxcurr = self.integrators[l2]._idxcurr

        l1sys, l2sys = self.integrators[l1].system, self.integrators[l2].system

        # Prevsoln is used as temporal storage at l1
        rtemp = 0 if l1idxcurr == 1 else 1

        # rtemp = R = -∇·f - dQ/dt
        self.integrator.system.rhs(self.tcurr, l1idxcurr, rtemp)

        # rtemp = -d = R - r at lower levels
        if l1 != self._order:
            self.integrator._add(1, rtemp, -1, self._mg_regidx[0])

        # Activate l2 system and get l2 regidx
        self.level = l2
        mg0, mg1 = self._mg_regidx

        # Restrict Q
        l1sys.eles_scal_upts_inb.active = l1idxcurr
        l2sys.eles_scal_upts_inb.active = l2idxcurr
        self.integrator._queue % self.mgproject(l1, l2)()

        # Restrict d and store to mg1
        l1sys.eles_scal_upts_inb.active = rtemp
        l2sys.eles_scal_upts_inb.active = mg1
        self.integrator._queue % self.mgproject(l1, l2)()

        # mg0 = R = -∇·f - dQ/dt
        self.integrator.system.rhs(self.tcurr, l2idxcurr, self._mg_regidx[0])

        # Compute the target residual r
        # mg0 = r = R + d
        self.integrator._add(1, self._mg_regidx[0], -1, self._mg_regidx[1])

        # Need to store the non-smoothed solution Q^ns for the correction
        # mg1 = Q^ns
        self.integrator._add(0, mg1, 1, l2idxcurr)

        # Restrict the physical stepper terms
        for i in range(self.integrator._stepper_nregs):
            l1sys.eles_scal_upts_inb.active = (
                self.integrators[l1]._stepper_regidx[i]
            )
            l2sys.eles_scal_upts_inb.active = (
                self.integrators[l2]._stepper_regidx[i]
            )
            self.integrator._queue % self.mgproject(l1, l2)()

    def prolongate(self, l1, l2):
        l1idxcurr  = self.integrators[l1]._idxcurr
        l2idxcurr = self.integrators[l2]._idxcurr

        l1sys, l2sys = self.integrators[l1].system, self.integrators[l2].system

        # Prevsoln is used as temporal storage at l2
        rtemp = 0 if l2idxcurr == 1 else 1

        # Correction with respect to the non-smoothed value from down-cycle
        # mg1 = Delta = Q^s - Q^ns
        self.integrator._add(-1, self._mg_regidx[1], 1, l1idxcurr)

        # Prolongate the correction and store to rtemp
        l1sys.eles_scal_upts_inb.active = self._mg_regidx[1]
        l2sys.eles_scal_upts_inb.active = rtemp
        self.integrator._queue % self.mgproject(l1, l2)()

        # Add the correction to the end quantity at l2
        # Q^m+1  = Q^s + Delta
        self.level = l2
        self.integrator._add(1, l2idxcurr, 1, rtemp)

    @property
    def _mg_regidx(self):
        if self.level == self._order:
            raise AttributeError('_mg_regidx not defined when'
                                 ' self.level == self._order')

        return self.integrator._aux_regidx[-2:]

    def pseudo_advance(self, tcurr, tout, dt):
        # Multigrid levels and step counts
        cycle, csteps = self.cycle, self.csteps

        self.tcurr = tcurr

        for i in range(self._maxniters):
            for l, m, n in it.zip_longest(cycle, cycle[1:], csteps):
                self.level = l
                # Set the number of smoothing steps at each level
                self.integrator.maxniters = self.integrator.minniters = n

                self.integrator.pseudo_advance(tcurr, tout, dt)

                if m is not None and l > m:
                    self.restrict(l, m)
                elif m is not None and l < m:
                    self.prolongate(l, m)

            # Convergence monitoring
            if self.mg_conv_mon(self.integrator, i, self._minniters):
                 break

        # Update the dual-time stepping banks
        self.finalise_mg_advance(self.integrator._idxcurr)
