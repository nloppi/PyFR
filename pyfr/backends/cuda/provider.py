# -*- coding: utf-8 -*-

import pkgutil

from mako.template import Template

import pycuda.compiler as compiler

from pyfr.util import memoize

class CudaKernelProvider(object):
    KERNEL_EXT = '.cu.mak'

    def _get_2d_grid_block(self, function, nrow, ncol):
        # TODO: Write a totally bitchin' method which uses info from the
        #       function to help compute an optimal block size
        block = (min(16, nrow), min(16, ncol), 1)
        grid = self._get_grid_for_block(block, nrow, ncol)
        return grid, block

    def _get_grid_for_block(self, block, nrow, ncol=1):
        return (nrow + (-nrow % block[0])) // block[0],\
               (ncol + (-ncol % block[1])) // block[1]

    @memoize
    def _get_module(self, module, tplparams={}, nvccopts=None):
        # Get the template file
        tpl = pkgutil.get_data(__name__, 'kernels/' + module + self.KERNEL_EXT)

        # Render the template
        mod = Template(tpl).render(**tplparams)

        # Compile
        return compiler.SourceModule(mod, options=nvccopts)

    @memoize
    def _get_function(self, module, function, argtypes, tplparams={},
                      nvccopts=None):
        # Compile/retrieve the module
        mod = self._get_module(module, tplparams, nvccopts)

        # Get a reference to the function
        func = mod.get_function(function)

        # Prepare it for execution
        return func.prepare(argtypes)
