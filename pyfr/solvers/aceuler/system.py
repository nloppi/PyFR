# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.aceuler.elements import ACEulerElements
from pyfr.solvers.aceuler.inters import (ACEulerIntInters, ACEulerMPIInters,
                                       ACEulerBaseBCInters)


class ACEulerSystem(BaseAdvectionSystem):
    name = 'ac-euler'

    elementscls = ACEulerElements
    intinterscls = ACEulerIntInters
    mpiinterscls = ACEulerMPIInters
    bbcinterscls = ACEulerBaseBCInters
