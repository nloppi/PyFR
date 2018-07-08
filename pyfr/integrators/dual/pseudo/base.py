# -*- coding: utf-8 -*-

from collections import defaultdict

from pyfr.integrators.basecommon import BaseCommon
from pyfr.util import proxylist

class BaseDualPseudoIntegrator(BaseCommon, object):
    formulation = 'dual'
    aux_nregs = 0

    def __init__(self, backend, systemcls, rallocs, mesh,
                 initsoln, cfg, tcoeffs):
        # Sanity checks
        if self._controller_needs_lerrest and not self._stepper_has_lerrest:
            raise TypeError('Incompatible stepper/controller combination')

        self.backend = backend
        self.rallocs = rallocs
        self.isrestart = initsoln is not None
        self.cfg = cfg

        sect = 'solver-time-integrator'

        self._dtaumin = 1.0e-12
        self._dtau = self.cfg.getfloat(sect, 'pseudo-dt')

        self._maxniters = self.cfg.getint(sect, 'pseudo-niters-max', 0)
        self._minniters = self.cfg.getint(sect, 'pseudo-niters-min', 0)

        if self._maxniters < self._minniters:
            raise ValueError('The maximum number of pseudo-iterations must '
                             'be greater than or equal to the minimum')

        self._pseudo_residtol= self.cfg.getfloat(sect, 'pseudo-resid-tol')
        self._pseudo_norm = self.cfg.get(sect, 'pseudo-resid-norm', 'l2')

        if self._pseudo_norm not in {'l2', 'uniform'}:
            raise ValueError('Invalid pseudo-residual norm')

        # Amount of temp storage required by physical stepper
        self._phys_stepper_nregs = len(tcoeffs) - 1

        # Determine the amount of temp storage required in total
        self.nreg = (self._pseudo_stepper_nregs + self._phys_stepper_nregs
                     + self.aux_nregs)

        # Check if multip integrator
        self.multip = 'multip' in self.name

        # BDF coefficients
        self._dual_time_source = tcoeffs

        # Construct the relevant mesh partition
        self.system = systemcls(backend, rallocs, mesh, initsoln,
                                nreg=self.nreg, cfg=cfg)

        # Storage for register banks and current index
        self._init_reg_banks(backend, self.system)

        # Additional registers for multip if present
        self._multip_regidx = self._regidx[-2:] if self.aux_nregs is not 0 else None

        # Get a queue for subclasses to use
        self._queue = backend.queue()

        # Global degree of freedom count
        self._gndofs = self._get_gndofs(self.system)

        elementscls = self.system.elementscls
        self._subdims = [elementscls.convarmap[self.system.ndims].index(v)
                         for v in elementscls.dualcoeffs[self.system.ndims]]

        # Pointwise kernels for integrator
        self.ikernels = {}
        self.intgkernels = defaultdict(proxylist)

    @property
    def _stepper_regidx(self):
        return self._regidx[:self._pseudo_stepper_nregs]

    @property
    def _source_regidx(self):
        pnreg, dsrc = self._pseudo_stepper_nregs, self._dual_time_source
        return self._regidx[pnreg:pnreg + self._phys_stepper_nregs]

    def _dual_time_source(self):
        pass

    def finalise_pseudo_advance(self, currsoln):
        pnreg = self._pseudo_stepper_nregs
        dtsnreg = len(self._dual_time_source)

        # Rotate the source registers to the right by one
        self._regidx[pnreg:pnreg + dtsnreg - 1] = (
            self._source_regidx[-1:] + self._source_regidx[:-1]
        )

        # Copy the current soln into the first source register
        self._add(0, self._regidx[pnreg], 1, currsoln)
