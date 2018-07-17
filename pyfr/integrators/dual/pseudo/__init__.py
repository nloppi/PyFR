# -*- coding: utf-8 -*-

import re

from pyfr.integrators.dual.pseudo.pseudocontrollers import BaseDualPseudoController
from pyfr.integrators.dual.pseudo.pseudosteppers import BaseDualPseudoStepper
from pyfr.integrators.dual.pseudo.multip import DualMultiPIntegrator
from pyfr.util import subclass_where


def get_pseudo_integrator(backend, systemcls, rallocs, mesh,
                          initsoln, cfg, tcoeffs):
    if 'solver-dual-time-integrator-multip' in cfg.sections():
        sect = 'solver-dual-time-integrator-multip'

        # Get the solver order
        order = cfg.getint('solver', 'order')

        # Get the multigrid cycle
        cycle, csteps = zip(*cfg.getliteral(sect, 'cycle'))
        levels = sorted(set(cycle), reverse=True)

        # Highest multip level
        pn = cfg.get('solver-time-integrator', 'pseudo-scheme')
        pc = subclass_where(BaseDualPseudoStepper, pseudo_stepper_name=pn)
        cn = cfg.get('solver-time-integrator', 'pseudo-controller')
        cc = subclass_where(BaseDualPseudoController, pseudo_controller_name=cn)
        cc_none = subclass_where(BaseDualPseudoController,
                                 pseudo_controller_name='none')
        # Get multip types
        multip_types = {}
        for l in levels:
            # No pseudo-control at lower levels
            cc = cc if l == order else cc_none
            bases = [(cn, cc), (pn, pc)]
            name = '_'.join(['dualP' + str(l)] + list(bn for bn, bc in bases) +
                            ['pseudointegrator'])
            multip_types[l] = type(name, tuple(bc for bn, bc in bases),
                                   dict(name=name))

        return DualMultiPIntegrator(backend, systemcls, rallocs, mesh,
                                    initsoln, cfg, tcoeffs, multip_types)
    else:
        cn = cfg.get('solver-time-integrator', 'pseudo-controller')
        pn = cfg.get('solver-time-integrator', 'pseudo-scheme')

        cc = subclass_where(BaseDualPseudoController, pseudo_controller_name=cn)
        pc = subclass_where(BaseDualPseudoStepper, pseudo_stepper_name=pn)

        bases = [(cn, cc), (pn, pc)]

        # Determine the integrator name
        name = '_'.join(['dual'] + list(bn for bn, bc in bases)
                        + ['pseudointegrator'])
        name = re.sub('(?:^|_|-)([a-z])', lambda m: m.group(1).upper(), name)

        pseudointegrator = type(name, tuple(bc for bn, bc in bases),
                                dict(name=name))

        # Construct and return an instance of this new integrator class
        return pseudointegrator(backend, systemcls, rallocs, mesh,
                                initsoln, cfg, tcoeffs)
