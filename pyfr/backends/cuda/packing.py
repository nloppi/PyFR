# -*- coding: utf-8 -*-

import pycuda.driver as cuda
from pycuda.compiler import SourceModule

from pyfr.backends.cuda.provider import CudaKernelProvider
from pyfr.backends.cuda.queue import CudaComputeKernel, CudaMPIKernel
from pyfr.backends.cuda.types import CudaMPIMatrix, CudaMPIView

from pyfr.util import npdtype_to_ctype, npdtype_to_mpitype

class CudaPackingKernels(CudaKernelProvider):
    def __init__(self, backend):
        pass

    def _packmodopts(self, mpiview):
        return dict(dtype=npdtype_to_ctype(mpiview.mpimat.dtype),
                    vlen=mpiview.view.vlen)

    def _packunpack_mpimat(self, op, mpimat):
        class PackUnpackKernel(CudaComputeKernel):
            if op == 'pack':
                def __call__(self, scomp, scopy):
                    cuda.memcpy_dtoh_async(mpimat.hdata, mpimat.data, scomp)
            else:
                def __call__(self, scomp, scopy):
                    cuda.memcpy_htod_async(mpimat.data, mpimat.hdata, scomp)

        return PackUnpackKernel()

    def _packunpack_mpiview(self, op, mpiview):
        # An MPI view is simply a regular view plus an MPI matrix
        v, m = mpiview.view, mpiview.mpimat

        # Get the CUDA pack/unpack kernel from the pack module
        fn = self._get_function('pack', op + '_view', 'iiPPPiii',
                                tplparams=self._packmodopts(mpiview))

        # We need all of the L1 we can get when packing
        fn.set_cache_config(cuda.func_cache.PREFER_L1)

        # Compute the grid and thread-block size
        grid, block = self._get_2d_grid_block(fn, v.nrow, v.ncol)

        # Create a CUDA event
        event = cuda.Event(cuda.event_flags.DISABLE_TIMING)

        class ViewPackUnpackKernel(CudaComputeKernel):
            def __call__(self, scomp, scopy):
                # If we are unpacking then copy the host buffer to the GPU
                if op == 'unpack':
                    cuda.memcpy_htod_async(m.data, m.hdata, scopy)
                    event.record(scopy)
                    scomp.wait_for_event(event)

                # Call the CUDA kernel (pack or unpack)
                fn.prepared_async_call(grid, block, scomp, v.nrow, v.ncol,
                                       v.mapping.data, v.strides.data, m.data,
                                       v.mapping.leaddim, v.strides.leaddim,
                                       m.leaddim)

                # If we have been packing then copy the GPU buffer to the host
                if op == 'pack':
                    event.record(scomp)
                    scopy.wait_for_event(event)
                    cuda.memcpy_dtoh_async(m.hdata, m.data, scopy)

        return ViewPackUnpackKernel()

    def _packunpack(self, op, mv):
        if isinstance(mv, CudaMPIMatrix):
            return self._packunpack_mpimat(op, mv)
        elif isinstance(mv, CudaMPIView):
            return self._packunpack_mpiview(op, mv)
        else:
            raise TypeError('Can only pack MPI views and MPI matrices')

    def _sendrecv(self, mv, mpipreqfn, pid, tag):
        # If we are an MPI view then extract the MPI matrix
        mpimat = mv.mpimat if isinstance(mv, CudaMPIView) else mv

        # Determine the MPI data type of the matrix
        mpitype = npdtype_to_mpitype(mpimat.dtype)

        # Create a persistent MPI request to send/recv the pack
        preq = mpipreqfn((mpimat.hdata, mpitype), pid, tag)

        class SendRecvPackKernel(CudaMPIKernel):
            def __call__(self, reqlist):
                # Start the request and append us to the list of requests
                preq.Start()
                reqlist.append(preq)

        return SendRecvPackKernel()

    def pack(self, mv):
        return self._packunpack('pack', mv)

    def send_pack(self, mv, mpicomm, pid, tag):
        return self._sendrecv(mv, mpicomm.Send_init, pid, tag)

    def recv_pack(self, mv, mpicomm, pid, tag):
        return self._sendrecv(mv, mpicomm.Recv_init, pid, tag)

    def unpack(self, mv):
        return self._packunpack('unpack', mv)
