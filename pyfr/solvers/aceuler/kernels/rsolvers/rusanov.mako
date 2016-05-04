# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.aceuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];

    ${pyfr.expand('inviscid_flux', 'ul', 'fl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr')};

% for i in range(ndims):
   vl[${i}] = ul[${i + 1}];
   vr[${i}] = ur[${i + 1}];
% endfor

    // Normal of the average interface velocity
    fpdtype_t nv = 0.5*${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Estimate the maximum wave speed
    fpdtype_t a_1 = fabs(nv - sqrt(nv*nv+${c['ac-zeta']}));
    fpdtype_t a_2 = fabs(nv + sqrt(nv*nv+${c['ac-zeta']}));
    fpdtype_t a_3 = fabs(nv);
    fpdtype_t a = max(max(a_1, a_2), a_3);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + 0.5*a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
