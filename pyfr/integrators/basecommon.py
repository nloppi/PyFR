# -*- coding: utf-8 -*-

from pyfr.util import proxylist, memoize
from pyfr.mpiutil import get_comm_rank_root, get_mpi

class BaseCommon(object):
    def _init_reg_banks(self, backend, system):
        self._regs, self._regidx = [], list(range(self.nreg))
        self._idxcurr = 0

        # Create a proxylist of matrix-banks for each storage register
        for i in self._regidx:
            self._regs.append(
                proxylist([backend.matrix_bank(em, i)
                           for em in system.ele_banks])
            )

    def _get_kernels(self, name, nargs, **kwargs):
        # Transpose from [nregs][neletypes] to [neletypes][nregs]
        transregs = zip(*self._regs)

        # Generate an kernel for each element type
        kerns = proxylist([])
        for tr in transregs:
            kerns.append(self.backend.kernel(name, *tr[:nargs], **kwargs))

        return kerns

    def _prepare_reg_banks(self, *bidxes):
        for reg, ix in zip(self._regs, bidxes):
            reg.active = ix

    def _get_gndofs(self, system):
        comm, rank, root = get_comm_rank_root()

        # Get the number of degrees of freedom in this partition
        ndofs = sum(system.ele_ndofs)

        # Sum to get the global number over all partitions
        return comm.allreduce(ndofs, op=get_mpi('sum'))

    @memoize
    def _get_axnpby_kerns(self, n, subdims=None):
        return self._get_kernels('axnpby', nargs=n, subdims=subdims)

    def _add(self, *args):
        # Get a suitable set of axnpby kernels
        axnpby = self._get_axnpby_kerns(len(args) // 2)

        # Bank indices are in odd-numbered arguments
        self._prepare_reg_banks(*args[1::2])

        # Bind and run the axnpby kernels
        self._queue % axnpby(*args[::2])
