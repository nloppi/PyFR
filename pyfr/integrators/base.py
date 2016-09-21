# -*- coding: utf-8 -*-

from abc import ABCMeta, abstractmethod, abstractproperty
from collections import deque
import re
import time

import numpy as np

from pyfr.inifile import Inifile
from pyfr.mpiutil import get_comm_rank_root, get_mpi
from pyfr.plugins import get_plugin
from pyfr.util import memoize, proxylist


class BaseIntegrator(object, metaclass=ABCMeta):
    def __init__(self, backend, systemcls, rallocs, mesh, initsoln, cfg):
        self.backend = backend
        self.rallocs = rallocs
        self.isrestart = initsoln is not None
        self.cfg = cfg
        self.prevcfgs = {f: initsoln[f] for f in initsoln or []
                         if f.startswith('config-')}

        # Ensure the system is compatible with our formulation
        if self.formulation not in systemcls.elementscls.formulations:
            raise RuntimeError(
                'System {0} does not support time stepping formulation {1}'
                .format(systemcls.name, self.formulation)
            )

        # Start time
        self.tstart = cfg.getfloat('solver-time-integrator', 'tstart', 0.0)
        self.tend = cfg.getfloat('solver-time-integrator', 'tend')

        # Current time; defaults to tstart unless restarting
        if self.isrestart:
            stats = Inifile(initsoln['stats'])
            self.tcurr = stats.getfloat('solver-time-integrator', 'tcurr')
        else:
            self.tcurr = self.tstart

        # List of target times to advance to
        self.tlist = deque([self.tend])

        # Accepted and rejected step counters
        self.nacptsteps = 0
        self.nrjctsteps = 0
        self.nacptchain = 0

        # Current and minimum time steps
        self._dt = self.cfg.getfloat('solver-time-integrator', 'dt')
        self.dtmin = 1.0e-12

        # Determine the amount of temp storage required by this method
        nreg = self._stepper_nregs

        # These are now lists. If no multigrid they hold only one element
        self.system = []
        self._regs = []
        self._regidx = []
        self._idxcurr = []

        if self.cfg.get('solver-time-integrator', 'controller') == 'mg':
            # Multigrid requires 3 additional registers
            nreg += 3
            self.leveliters = self.cfg.getintlist('solver-time-integrator',
                                                  'cycle')
            self.levels = len(self.leveliters)

            # Initialise multiple systems for multigird
            for level in range(self.levels):
                self.system.append(systemcls(backend, rallocs, mesh, initsoln,
                                              nreg, cfg, level))
                self._regs.append(self._get_reg_banks(nreg, level)[0])
                self._regidx.append(self._get_reg_banks(nreg, level)[1])
                self._idxcurr.append(0)

        else:
            # Construct the relevant mesh partition
            self._idxcurr.append(0)
            self.system.append(systemcls(backend, rallocs, mesh,
                                         initsoln, nreg, cfg))
            # Storage register banks
            self._regs.append(self._get_reg_banks(nreg, level=0)[0])
            self._regidx.append(self._get_reg_banks(nreg, level=0)[1])

        # Extract the UUID of the mesh (to be saved with solutions)
        self.mesh_uuid = mesh['mesh_uuid']

        # Get a queue for subclasses to use
        self._queue = backend.queue()

        # Global degree of freedom count
        self._gndofs = self._get_gndofs()

        # Solution cache
        self._curr_soln = None

        # Add kernel cache
        self._axnpby_kerns = {}

        # Record the starting wall clock time
        self._wstart = time.time()

        # Event handlers for advance_to
        self.completed_step_handlers = proxylist(self._get_plugins())

        # Delete the memory-intensive elements map from the system
        #del self.system.ele_map

    def _get_reg_banks(self, nreg, level):
        regs, regidx = [], list(range(nreg))

        # Create a proxylist of matrix-banks for each storage register
        for i in regidx:
            regs.append(
                proxylist([self.backend.matrix_bank(em, i)
                           for em in self.system[level].ele_banks])
            )

        return regs, regidx

    def _get_gndofs(self):
        comm, rank, root = get_comm_rank_root()

        # Get the number of degrees of freedom in this partition
        ndofs = sum(self.system[0].ele_ndofs)

        # Sum to get the global number over all partitions
        return comm.allreduce(ndofs, op=get_mpi('sum'))

    def _get_plugins(self):
        plugins = []

        for s in self.cfg.sections():
            m = re.match('soln-plugin-(.+?)(?:-(.+))?$', s)
            if m:
                cfgsect, name, suffix = m.group(0), m.group(1), m.group(2)

                # Instantiate
                plugins.append(get_plugin(name, self, cfgsect, suffix))

        return plugins

    def _get_kernels(self, name, nargs, level=0, **kwargs):
        # Transpose from [nregs][neletypes] to [neletypes][nregs]
        transregs = zip(*self._regs[0])

        # Generate an kernel for each element type
        kerns = proxylist([])
        for tr in transregs:
            kerns.append(self.backend.kernel(name, *tr[:nargs], **kwargs))

        return kerns

    def _prepare_reg_banks(self, *bidxes, level=0):
        for reg, ix in zip(self._regs[0], bidxes):
            reg.active = ix

    @memoize
    def _get_axnpby_kerns(self, n, level=0, subdims=None):
        return self._get_kernels('axnpby', nargs=n, level=level, subdims=subdims)

    def _add(self, *args, level=0):
        # Get a suitable set of axnpby kernels
        axnpby = self._get_axnpby_kerns(len(args) // 2, level=level)

        # Bank indices are in odd-numbered arguments
        self._prepare_reg_banks(*args[1::2], level=level)

        # Bind and run the axnpby kernels
        self._queue % axnpby(*args[::2])

    def call_plugin_dt(self, dt):
        ta = self.tlist
        tb = deque(np.arange(self.tcurr, self.tend, dt).tolist())

        self.tlist = tlist = deque()

        # Merge the current and new time lists
        while ta and tb:
            t = ta.popleft() if ta[0] < tb[0] else tb.popleft()
            if not tlist or t - tlist[-1] > self.dtmin:
                tlist.append(t)

        tlist.extend(ta)
        tlist.extend(tb)

    @property
    def soln(self):
        # If we do not have the solution cached then fetch it
        if not self._curr_soln:
            self._curr_soln = self.system[0].ele_scal_upts(self._idxcurr[0])

        return self._curr_soln

    @abstractmethod
    def step(self, t, dt):
        pass

    @abstractmethod
    def advance_to(self, t):
        pass

    @abstractproperty
    def _stepper_nfevals(self):
        pass

    @abstractproperty
    def _stepper_nregs(self):
        pass

    @abstractproperty
    def _stepper_order(self):
        pass

    def run(self):
        for t in self.tlist:
            self.advance_to(t)

    @property
    def nsteps(self):
        return self.nacptsteps + self.nrjctsteps

    def collect_stats(self, stats):
        wtime = time.time() - self._wstart

        # Rank allocation
        stats.set('backend', 'rank-allocation',
                  ','.join(str(r) for r in self.rallocs.mprankmap))

        # Simulation and wall clock times
        stats.set('solver-time-integrator', 'tcurr', self.tcurr)
        stats.set('solver-time-integrator', 'wall-time', wtime)

        # Step counts
        stats.set('solver-time-integrator', 'nsteps', self.nsteps)
        stats.set('solver-time-integrator', 'nacptsteps', self.nacptsteps)
        stats.set('solver-time-integrator', 'nrjctsteps', self.nrjctsteps)

    @property
    def cfgmeta(self):
        cfg = self.cfg.tostr()

        if self.prevcfgs:
            ret = dict(self.prevcfgs, config=cfg)

            if cfg != ret['config-' + str(len(self.prevcfgs) - 1)]:
                ret['config-' + str(len(self.prevcfgs))] = cfg

            return ret
        else:
            return {'config': cfg, 'config-0': cfg}

    @memoize
    def prolrest(self, l1, l2, source=None):
        prolrestkerns = proxylist([])

        for etype in self.system[0].ele_types:
            l1sys = self.system[l1]
            l2sys = self.system[l2]

            l1eles = l1sys.ele_map[etype]
            l2eles = l2sys.ele_map[etype]

            # fix
            if l2 < l1:
                prolrestkerns.append(
                    self.backend.kernel('mul',
                                        self._get_prolong(l1eles, l2eles),
                                        l1eles.scal_upts_inb,
                                        out=l2eles.scal_upts_inb)
                )

            elif l2 > l1:
                prolrestkerns.append(
                    self.backend.kernel('mul',
                                        self._get_prolong(l1eles, l2eles), #NOW UNDERSAMPLING INSTEAD OF VDM CHOP
                                        l1eles.scal_upts_inb,
                                        out=l2eles.scal_upts_inb)
                )

            elif source:
                prolrestkerns.append(
                    self.backend.kernel('mul',
                                        self._get_prolong(l1eles, l2eles),
                                        l1eles.scal_upts_inb,
                                        out=l2eles.scal_upts_inb)
                )
            else:
                print('Restriction/Prolongation must be done '
                  ' between two different polynomial levels')

        return prolrestkerns

    def restrict(self, l1, l2):
        rhs1 = self.system[l1].rhs
        rhs2 = self.system[l2].rhs
        idxcurr = self._idxcurr
        regidx = self._regidx

        mgregidx = self._MG_regidx
        add_with_dts = self._add_with_dts

        # -R temporarility to last mgstepper index of l1
        rhs1(self.tcurr, idxcurr[l1], regidx[l1][mgregidx[1]])
        if l1 == 0:
            l = 0
        else:
            l = 1

        # d = r - (rhs - source(r1)) to the second MG index of l1. check signs
        add_with_dts(-1, regidx[l1][mgregidx[1]], l, regidx[l1][mgregidx[0]],
                     c=-1/self._dt,  level=l1)

        # Restrict Q and d, store Q^start to the last stepper index of l2,
        # compute r
        self.system[l1].eles_scal_upts_inb.active = idxcurr[l1]
        self.system[l2].eles_scal_upts_inb.active = idxcurr[l2]
        self._queue % self.prolrest(l1, l2)()

        self._add(0.0, regidx[l2][mgregidx[-1]], 1, idxcurr[l2], level=l2)

        self.system[l1].eles_scal_upts_inb.active = regidx[l1][mgregidx[1]]
        self.system[l2].eles_scal_upts_inb.active = regidx[l2][mgregidx[1]]
        self._queue % self.prolrest(l1, l2)()

        # -R(Q^s) temporarility to MG regidx0
        rhs2(0.0, idxcurr[l2], regidx[l2][mgregidx[0]])

        # r = -RHS - dtssrc + d
        add_with_dts(1, regidx[l2][mgregidx[0]], 1, regidx[l2][mgregidx[1]],
                     c=1/self._dt, level=l2)

        # Source terms at new solnpts
        for srcidx in self._source_regidx:
            self.system[l1].eles_scal_upts_inb.active = self._regidx[l1][srcidx]
            self.system[l2].eles_scal_upts_inb.active = self._regidx[l2][srcidx]
            self.prolrest(l1, l2, source=True)
            self._queue % self.prolrest(l1, l2)()

    def prolongate(self, l1, l2):
        idxcurrs = self._idxcurr
        regidx = self._regidx
        mgregidx = self._MG_regidx

        # Calculate the correction
        self._add(0, regidx[l1][mgregidx[1]], 1, idxcurrs[l1], -1,
                  regidx[l1][mgregidx[-1]], level=l1)

        # Prolongation for the correction
        self.system[l1].eles_scal_upts_inb.active = regidx[l1][mgregidx[1]]
        self.system[l2].eles_scal_upts_inb.active = regidx[l2][mgregidx[1]]
        self._queue % self.prolrest(l1, l2)()

        # Add the correction to the end quantity at l2
        self._add(1, idxcurrs[l2], 1, regidx[l2][mgregidx[1]], level=l2)

    def _get_restrict(self, name, l1eles, l2eles):
        l1b = l1eles.basis
        l2b = l2eles.basis
        ord = l1b.order
        ub1 = l1b.ubasis
        ub2 = l2b.ubasis

        if name == 'quad':
            a = np.eye(l2b.nupts, l1b.nupts)
            for i in range(1, ord):
                a[(ord*i):ord*(i + 1)] = np.roll(a[(ord*i):ord*(i + 1)], i,
                                                 axis=1)

        if name == 'tri':
            a = np.eye(ub2.nupts, ub1.nupts)
            j = 0
            for i in range(ord - 1):
                j += ord - i
                a[j:j+(ord - i - 1)] = np.roll(a[j:j + (ord - i - 1)], i + 1,
                                               axis=1)

        mat = np.dot(ub2.vdm, np.dot(a, np.linalg.inv(ub1.vdm)))
        return self.backend.const_matrix(mat, tags={'align'})

    def _get_prolong(self, l1eles, l2eles):
        ub1 = l1eles.basis.ubasis
        l2b = l2eles.basis

        mat = ub1.nodal_basis_at(l2b.upts)
        return self.backend.const_matrix(mat, tags={'align'})
