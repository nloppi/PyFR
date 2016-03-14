.. highlightlang:: python

***************
Developer Guide
***************

======================================
A Brief Overview of the PyFR Framework
======================================

Where to Start
--------------

The symbolic link :code:`pyfr.scripts.pyfr` points to the script
:code:`pyfr.scripts.main`, which is where it all starts! Specifically,
the function :code:`process_run` calls the function
:code:`_process_common`, which in turn calls the function
:code:`get_solver`, returning an Integrator -- a composite of a
`Controller`_ and a `Stepper`_. The Integrator has a method named
:code:`run`, which is then called to run the simulation.

Controller
----------

A `Controller`_ acts to advance the simulation in time. Specifically, a
`Controller`_ has a method named :code:`advance_to` which advances a
`System`_ to a specified time. There are two types of `Controller`_
available in PyFR |release|:

.. autoclass:: pyfr.integrators.controllers.NoneController
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.integrators.controllers.PIController
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

Types of `Controller`_ are related via the following inheritance
diagram:

.. inheritance-diagram:: pyfr.integrators.controllers
    :parts: 1

Stepper
-------

A `Stepper`_ acts to advance the simulation by a single time-step.
Specifically, a `Stepper`_ has a method named :code:`step` which
advances a `System`_ by a single time-step. There are five types of
`Stepper`_ available in PyFR |release|:

.. autoclass:: pyfr.integrators.steppers.EulerStepper
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.integrators.steppers.RK4Stepper
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.integrators.steppers.RK34Stepper
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.integrators.steppers.RK45Stepper
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.integrators.steppers.TVDRK3Stepper
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

Types of `Stepper`_ are related via the following inheritance diagram:

.. inheritance-diagram:: pyfr.integrators.steppers
    :parts: 1

System
------

A `System`_ holds information/data for the system, including
`Elements`_, `Interfaces`_, and the `Backend`_ with which the
simulation is to run. A `System`_ has a method named :code:`rhs`, which
obtains the divergence of the flux (the 'right-hand-side') at each
solution point. The method :code:`rhs` invokes various kernels which
have been pre-generated and loaded into queues. A `System`_ also has a
method named :code:`_gen_kernels` which acts to generate all the
kernels required by a particular `System`_. A kernel is an instance of
a 'one-off' class with a method named :code:`run` that implements the
required kernel functionality. Individual kernels are produced by a
kernel provider. PyFR |release| has various types of kernel provider. A
`Pointwise Kernel Provider`_ produces point-wise kernels such as
Riemann solvers and flux functions etc. These point-wise kernels are
specified using an in-built platform-independent templating language
derived from `Mako <http://www.makotemplates.org/>`_, henceforth
referred to as `PyFR-Mako`_. There are two types of `System`_ available
in PyFR |release|:

.. autoclass:: pyfr.solvers.euler.system.EulerSystem
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.solvers.navstokes.system.NavierStokesSystem
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

Types of `System`_ are related via the following inheritance diagram:

.. inheritance-diagram:: pyfr.solvers.navstokes.system
                         pyfr.solvers.euler.system
    :parts: 1

Elements
--------

An `Elements`_ holds information/data for a group of elements. There are
two types of `Elements`_ available in PyFR |release|:

.. autoclass:: pyfr.solvers.euler.elements.EulerElements
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.solvers.navstokes.elements.NavierStokesElements
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

Types of `Elements`_ are related via the following inheritance diagram:

.. inheritance-diagram:: pyfr.solvers.navstokes.elements
                         pyfr.solvers.euler.elements
    :parts: 1

Interfaces
----------

An `Interfaces`_ holds information/data for a group of interfaces. There
are four types of (non-boundary) `Interfaces`_ available in PyFR
|release|:

.. autoclass:: pyfr.solvers.euler.inters.EulerIntInters
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.solvers.euler.inters.EulerMPIInters
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.solvers.navstokes.inters.NavierStokesIntInters
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.solvers.navstokes.inters.NavierStokesMPIInters
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

Types of (non-boundary) `Interfaces`_ are related via the following
inheritance diagram:

.. inheritance-diagram:: pyfr.solvers.navstokes.inters.NavierStokesMPIInters
                         pyfr.solvers.navstokes.inters.NavierStokesIntInters
                         pyfr.solvers.euler.inters.EulerMPIInters
                         pyfr.solvers.euler.inters.EulerIntInters
    :parts: 1

Backend
-------

A `Backend`_ holds information/data for a backend. There are three types
of `Backend`_ available in PyFR |release|:

.. autoclass:: pyfr.backends.cuda.base.CUDABackend
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.opencl.base.OpenCLBackend
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.openmp.base.OpenMPBackend
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

Types of `Backend`_ are related via the following inheritance diagram:

.. inheritance-diagram:: pyfr.backends.cuda.base
                         pyfr.backends.opencl.base
                         pyfr.backends.openmp.base
    :parts: 1

Pointwise Kernel Provider
-------------------------

A `Pointwise Kernel Provider`_ produces point-wise kernels.
Specifically, a `Pointwise Kernel Provider`_ has a method named
:code:`register`, which adds a new method to an instance of a
`Pointwise Kernel Provider`_. This new method, when called, returns a
kernel. A kernel is an instance of a 'one-off' class with a method
named :code:`run` that implements the required kernel functionality.
The kernel functionality itself is specified using `PyFR-Mako`_. Hence,
a `Pointwise Kernel Provider`_ also has a method named
:code:`_render_kernel`, which renders `PyFR-Mako`_ into low-level
platform-specific code. The :code:`_render_kernel` method first sets
the context for Mako (i.e. details about the `Backend`_ etc.) and then
uses Mako to begin rendering the `PyFR-Mako`_ specification. When Mako
encounters a :code:`pyfr:kernel` an instance of a `Kernel Generator`_
is created, which is used to render the body of the
:code:`pyfr:kernel`. There are three types of `Pointwise Kernel
Provider`_ available in PyFR |release|:

.. autoclass:: pyfr.backends.cuda.provider.CUDAPointwiseKernelProvider
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.opencl.provider.OpenCLPointwiseKernelProvider
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.openmp.provider.OpenMPPointwiseKernelProvider
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

Types of `Pointwise Kernel Provider`_ are related via the following
inheritance diagram:

.. inheritance-diagram:: pyfr.backends.openmp.provider
                         pyfr.backends.cuda.provider
                         pyfr.backends.opencl.provider
                         pyfr.backends.base.kernels.BasePointwiseKernelProvider
    :parts: 1

Kernel Generator
----------------

A `Kernel Generator`_ renders the `PyFR-Mako`_ in a :code:`pyfr:kernel`
into low-level platform-specific code. Specifically, a `Kernel
Generator`_ has a method named :code:`render`, which applies `Backend`_
specific regex and adds `Backend`_ specific 'boiler plate' code to
produce the low-level platform-specific source -- which is compiled,
linked, and loaded. There are three types of `Kernel Generator`_
available in PyFR |release|:

.. autoclass:: pyfr.backends.cuda.generator.CUDAKernelGenerator
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.opencl.generator.OpenCLKernelGenerator
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.openmp.generator.OpenMPKernelGenerator
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

Types of `Kernel Generator`_ are related via the following inheritance diagram:

.. inheritance-diagram:: pyfr.backends.cuda.generator.CUDAKernelGenerator
                         pyfr.backends.opencl.generator.OpenCLKernelGenerator
                         pyfr.backends.openmp.generator.OpenMPKernelGenerator
    :parts: 1

=========
Data Stucture
=========

Matrices
--------

Matrices in PyFR are handled with a backend specific
:class:`~pyfr.backends.base.types.Matrix` objects. The memory allocation for
matrices is done by `System`_ that gathers the necessary information to shape
the matrix, such as number of elements, solution points and variables.
Regardless of the dimensionality of the matrix,
:class:`~pyfr.backends.base.types.Matrix` class treats them
as two-dimensional, using a representation illustrated below for a
4-dimensional matrix with a shape
:code:`self.shape = (ndims=2, nupts, nvars=3, neles)`.

..  tikz:: [scale=0.6]
  \draw[shift={(0.5,0.5)}] (0,0) grid (19,10);
  \foreach \x in {6} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=gray] at (\x,\y) {};;}}
  \foreach \x in {12} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=gray] at (\x,\y) {};;}}
  \foreach \x in {18} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=gray] at (\x,\y) {};;}}
  \foreach \x in {19} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=green] at (\x,\y) {};;}}
  \node at (3,-0.15) {$\text{leadsubdim}$};
  \draw[<->] (0.5,-0.5) -- (6.5,-0.5);
  \node at (9,-0.95) {$\text{leaddim}$};
  \draw[<->] (0.5,-1.3) -- (19.5,-1.3);
  \draw[line width=1mm] (6.5, 0.5) -- (6.5, 10.5);
  \draw[line width=1mm] (12.5, 0.5) -- (12.5, 10.5);
  \draw[line width=1mm] (0.5, 5.5) -- (19.5, 5.5);
  \draw[dashed, line width=0.6mm, dash pattern=on 5pt off 3pt] (17.5, 0.5) -- (17.5, 10.5);
  \draw[-] (6,-0.1) -- (18,-0.1);
  \draw[->] (18,-0.1) -- (18,0.3);
  \draw[->] (6,-0.1) -- (6,0.3);
  \draw[->] (12,-0.1) -- (12,0.3);
  \node at (12,-0.5) {$\text{padding}$};
  \draw[<->] (0.5, 10.8) -- (5.5, 10.8);
  \node at (3,11.2) {$\text{neles}$};
  \draw[<->] (0.1, 10.5) -- (0.1, 5.5);
  \node at (-1,8) {$\text{nupts}$};
  \draw[->] (19,11.2) -- (19,10.7);
  \node at (18,11.5) {$\text{BLAS padding}$};

Padding is added to the subdimensions to adjust their lengths to be
divisible with the memory alignment requirement. Additional padding can
added at the end of the leaddimension to satisfy the alignment requirement of
some BLAS library methods.

Individual elements in the row major matrix are referenced with a formula
:math:`f(A,B,C,D) = A \cdot B \cdot leaddim + B \cdot leaddim + C \cdot leadsubdim + D`. Multiplication
kernels essentially loop over this formula using
instruction stream blocks/threads to calculate a group of elements in
parallel. The padding for the last column can be neglected in the
matrix operations as it exist only to fulfill the memory alignment
requirement and do not contain data that needs to be operated.

.. autoclass:: pyfr.backends.base.types.MatrixBase
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry


.. autoclass:: pyfr.backends.cuda.types.CUDAMatrixBase
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.opencl.types.OpenCLMatrixBase
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry


.. autoclass:: pyfr.backends.openmp.types.OpenMPMatrixBase
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.base.types.Matrix
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry


.. autoclass:: pyfr.backends.cuda.types.CUDAMatrix
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.opencl.types.OpenCLMatrix
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry


.. autoclass:: pyfr.backends.openmp.types.OpenMPMatrix
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. inheritance-diagram:: pyfr.backends.base.types.MatrixBase
                         pyfr.backends.cuda.types.CUDAMatrixBase
                         pyfr.backends.opencl.types.OpenCLMatrixBase
                         pyfr.backends.openmp.types.OpenMPMatrixBase
                         pyfr.backends.base.types.Matrix
                         pyfr.backends.cuda.types.CUDAMatrix
                         pyfr.backends.opencl.types.OpenCLMatrix
                         pyfr.backends.openmp.types.OpenMPMatrix
    :parts: 1



Matrix Banks
------------

:class:`~pyfr.backends.base.types.MatrixBank` objects are a combination of
several matrices having the same shape. Matrices in a
:class:`~pyfr.backends.base.types.MatrixBank` have an allocated index
that points to a register containing the Matrix. Individual
matrices in :class:`~pyfr.backends.base.types.MatrixBank` can be operated
by making them active, when the MatrixBank essentially represents
a single Matrix. A schematic figure below illustrates the
structure and usage of a :class:`~pyfr.backends.base.types.MatrixBank`

..  tikz:: [scale=1.0]
  \draw[shift={(0.5,0.5)}] (0,0) grid (4,1);
  \node at (1,1) {$[A]$};
  \node at (2,1) {$[B]$};
  \node at (3,1) {$[C]$};
  \node at (4,1) {$[D]$};
  \node at (0.1,2) {$regidx=[0$};
  \node at (2,2) {$1$};
  \node at (3,2) {$2$};
  \node at (4,2) {$3]$};


..  tikz:: [scale=1.0]
  \draw[shift={(0.5,0.5)}] (0,0) grid (4,1);
  \node[fill=green, minimum size=0.5cm] at (3,1) {};
  \node at (1,1) {$[A]$};
  \node at (2,1) {$[B]$};
  \node at (3,1) {$[C]$};
  \node at (4,1) {$[D]$};
  \node at (0.1,2) {\texttt{MatrixBank.active = 2}};
  \node at (1,0) {$\texttt{NewMatrix} = \texttt{MatrixBank} \times 3 = 3[C]$};

..  tikz:: [scale=1.0]
  \draw[shift={(0.5,0.5)}] (0,0) grid (4,1);
  \node[fill=green, minimum size=0.5cm] at (1,1) {};
  \node at (1,1) {$[A]$};
  \node at (2,1) {$[B]$};
  \node at (3,1) {$[C]$};
  \node at (4,1) {$[D]$};
  \node at (0.1,2) {\texttt{MatrixBank.active = 0}};
  \node at (2.3,0) {$\texttt{AnotherMatrix} = \texttt{MatrixBank} + \texttt{MatrixBank} = 2[A]$};


.. autoclass:: pyfr.backends.base.types.MatrixBank
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry


.. autoclass:: pyfr.backends.cuda.types.CUDAMatrixBank
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.opencl.types.OpenCLMatrixBank
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry


.. autoclass:: pyfr.backends.openmp.types.OpenMPMatrixBank
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. inheritance-diagram:: pyfr.backends.cuda.types.CUDAMatrixBank
                         pyfr.backends.opencl.types.OpenCLMatrixBank
                         pyfr.backends.openmp.types.OpenMPMatrixBank
    :parts: 1


Views
------------

PyFR uses Riemann solvers to calculate common interface for flux correction.
Information on both sides of the interface is needed and View class
is used to keep track and pair the numerical flux interface points
that correspond to each other.

..  tikz:: [scale=0.7]
  \draw (0,4) -- (4,4) -- (2,6) -- (0,4);
  \draw (0,0) -- (4,0) -- (4,4) -- (0,4) -- (0,0);
  \filldraw (2,3.85) circle (1.5pt);
  \filldraw (0.4,3.85) circle (1.5pt);
  \filldraw (3.6,3.85) circle (1.5pt);
  \filldraw (2,4.15) circle (1.5pt);
  \filldraw (0.4,4.15) circle (1.5pt);
  \filldraw (3.6,4.15) circle (1.5pt);
  \draw (1.8,3.7) -- (2.2,3.7) -- (2.2,4.3) -- (1.8,4.3) -- (1.8,3.7);
  \node at (2,3.3) {$lhs$};
  \node at (2,4.6) {$rhs$};

Let us consider two matrices containing the divergence of conservative variables
in all flux points for two element types. To generate a View matrices that holds only
the boundary values for left-hand side :math:`lhs` and right-hand side :math:`rhs`, first we need to get the
view permutation of the these nodes. The permutation
is arbitrary and results in an optimal memory access pattern for the :math:`LHS` of the interface.
To generate the View Matrix, we simply need to pass the permutation together with the
extent, for example::

            self._gen_perm(lhs, rhs)
            self._scal0_lhs = self._scal_view(lhs, 'get_scal_fpts_for_inter')
            self._scal0_rhs = self._scal_view(rhs, 'get_scal_fpts_for_inter')

The underlying methodology behind the extent, is passing several individual arguments
for the :class:`~pyfr.backends.base.types.View` class initialisation:

 * :code:`matmap` - A list of matrix ids(e.g. :code:`[trianglesmat_id, quadsmat_id]`) that we want
   to be viewed
 * :code:`rcmap`  - A list that contains row and column indices of the matrices
 * :code:`rcstride` (:code:`[[rstride], [cstrice]]`) - A list that defines how many
   matrix elements we need jump forward and downward to get to the next conservative
   variable and other dimension, respectively.
 * :code:`vshape` shape of the matrices for book keeping
 * :code:`tags` Matrix tags for book keeping


..  tikz:: [scale=0.4]
  \draw[shift={(0.5,0.5)}] (0,0) grid (19,10);
  \foreach \x in {6} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=gray] at (\x,\y) {};;}}
  \foreach \x in {12} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=gray] at (\x,\y) {};;}}
  \foreach \x in {18} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=gray] at (\x,\y) {};;}}
  \foreach \x in {19} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=green] at (\x,\y) {};;}}
  \draw[line width=1mm] (6.5, 0.5) -- (6.5, 10.5);
  \draw[line width=1mm] (12.5, 0.5) -- (12.5, 10.5);
  \draw[line width=1mm] (0.5, 5.5) -- (19.5, 5.5);
  \draw[dashed, line width=0.6mm, dash pattern=on 5pt off 3pt] (17.5, 0.5) -- (17.5, 10.5);
  \node[fill=white, minimum size=0.7cm] at (3,9.5) {};
  \node at (3,9.5) {$\nabla \rho$};
  \node[fill=white, minimum size=0.7cm] at (9,9.5) {};
  \node at (9,9.5) {$\nabla \rho u$};
  \node[fill=white, minimum size=0.7cm] at (15,9.5) {};
  \node at (15,9.5) {$\nabla E$};
  \node[right] at (23,5.5) {Matrix1: Triangles};
  \node[fill=white] at (9,7) {\texttt{cstride}};
  \node[fill=white] at (2.8,3) {\texttt{rstride}};
  \draw[->, line width=0.4mm ] (5,6) -- (11,6);
  \draw[->, line width=0.4mm ] (5,6) -- (5,1);


..  tikz:: [scale=0.4]
  \draw[shift={(0.5,0.5)}] (0,0) grid (19,10);
  \foreach \x in {6} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=gray] at (\x,\y) {};;}}
  \foreach \x in {12} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=gray] at (\x,\y) {};;}}
  \foreach \x in {18} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=gray] at (\x,\y) {};;}}
  \foreach \x in {19} {
  \foreach \y in {1,2,3,4,5,6,7,8,9,10} {
  \node[fill=green] at (\x,\y) {};;}}
  \draw[line width=1mm] (6.5, 0.5) -- (6.5, 10.5);
  \draw[line width=1mm] (12.5, 0.5) -- (12.5, 10.5);
  \draw[line width=1mm] (0.5, 5.5) -- (19.5, 5.5);
  \draw[dashed, line width=0.6mm, dash pattern=on 5pt off 3pt] (17.5, 0.5) -- (17.5, 10.5);
  \node[fill=white, minimum size=0.7cm] at (3,9.5) {};
  \node at (3,9.5) {$\nabla \rho$};
  \node[fill=white, minimum size=0.7cm] at (9,9.5) {};
  \node at (9,9.5) {$\nabla \rho u$};
  \node[fill=white, minimum size=0.7cm] at (15,9.5) {};
  \node at (15,9.5) {$\nabla E$};
  \node[right] at (23,5.5) {Matrix1: Quadrilaterals};


.. autoclass:: pyfr.backends.base.types.View
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.cuda.types.CUDAView
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry

.. autoclass:: pyfr.backends.opencl.types.OpenCLView
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry


.. autoclass:: pyfr.backends.openmp.types.OpenMPView
    :members:
    :undoc-members:
    :inherited-members:
    :private-members:
    :exclude-members: _abc_cache, _abc_negative_cache,
                      _abc_negative_cache_version, _abc_registry


.. inheritance-diagram:: pyfr.backends.cuda.types.CUDAView
                         pyfr.backends.opencl.types.OpenCLView
                         pyfr.backends.openmp.types.OpenMPView
    :parts: 1


=========
PyFR-Mako
=========

PyFR-Mako Kernels
-----------------

PyFR-Mako kernels are specifications of point-wise functionality that
can be invoked directly from within PyFR. They are opened with a header
of the form::

    <%pyfr:kernel name='kernel-name' ndim='data-dimensionality' [argument-name='argument-intent argument-attribute argument-data-type' ...]>

where

1. ``kernel-name`` --- name of kernel

    *string*

2. ``data-dimensionality`` --- dimensionality of data

    *int*

3. ``argument-name`` --- name of argument

    *string*

4. ``argument-intent`` --- intent of argument

    ``in`` | ``out`` | ``inout``

5. ``argument-attribute`` --- attribute of argument

    ``mpi`` | ``scalar`` | ``view``

6. ``argument-data-type`` --- data type of argument

    *string*

and are closed with a footer of the form::

     </%pyfr:kernel>

PyFR-Mako Macros
----------------

PyFR-Mako macros are specifications of point-wise functionality that
cannot be invoked directly from within PyFR, but can be embedded into
PyFR-Mako kernels. PyFR-Mako macros can be viewed as building blocks
for PyFR-mako kernels. They are opened with a header of the form::

    <%pyfr:macro name='macro-name' params='[parameter-name, ...]'>

where

1. ``macro-name`` --- name of macro

    *string*

2. ``parameter-name`` --- name of parameter

    *string*

and are closed with a footer of the form::

    </%pyfr:macro>

PyFR-Mako macros are embedded within a kernel using an expression of
the following form::

    ${pyfr.expand('macro-name', ['parameter-name', ...])};

where

1. ``macro-name`` --- name of the macro

    *string*

2. ``parameter-name`` --- name of parameter

    *string*

Syntax
------

Basic Functionality
^^^^^^^^^^^^^^^^^^^

Basic functionality can be expressed using a restricted subset of the C
programming language. Specifically, use of the following is allowed:

1. ``+,-,*,/`` --- basic arithmetic

2. ``sin, cos, tan`` --- basic trigonometric functions

3. ``exp`` --- exponential

4. ``pow`` --- power

5. ``fabs`` --- absolute value

6. ``output = ( condition ? satisfied : unsatisfied )`` --- ternary if

7. ``min`` --- minimum

8. ``max`` --- maximum

However, conditional if statements, as well as for/while loops, are
not allowed.

Expression Substitution
^^^^^^^^^^^^^^^^^^^^^^^

Mako expression substitution can be used to facilitate PyFR-Mako kernel
specification. A Python expression :code:`expression` prescribed thus
:code:`${expression}` is substituted for the result when the PyFR-Mako
kernel specification is interpreted at runtime.

Example::

        E = s[${ndims - 1}]

Conditionals
^^^^^^^^^^^^

Mako conditionals can be used to facilitate PyFR-Mako kernel
specification. Conditionals are opened with :code:`% if condition:` and
closed with :code:`% endif`. Note that such conditionals are evaluated
when the PyFR-Mako kernel specification is interpreted at runtime, they
are not embedded into the low-level kernel.

Example::

        % if ndims == 2:
            fout[0][1] += t_xx;     fout[1][1] += t_xy;
            fout[0][2] += t_xy;     fout[1][2] += t_yy;
            fout[0][3] += u*t_xx + v*t_xy + ${-c['mu']*c['gamma']/c['Pr']}*T_x;
            fout[1][3] += u*t_xy + v*t_yy + ${-c['mu']*c['gamma']/c['Pr']}*T_y;
        % endif

Loops
^^^^^

Mako loops can be used to facilitate PyFR-Mako kernel specification.
Loops are opened with :code:`% for condition:` and closed with :code:`%
endfor`. Note that such loops are unrolled when the PyFR-Mako kernel
specification is interpreted at runtime, they are not embedded into the
low-level kernel.

Example::

        % for i in range(ndims):
            rhov[${i}] = s[${i + 1}];
            v[${i}] = invrho*rhov[${i}];
        % endfor
