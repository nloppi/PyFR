# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, spngmuin, fout'>
    fpdtype_t divu = ${1.0/3.0}*(${' + '.join('grad_uin[{i}][{j}]'
                                              .format(i=i,j=i + 1)
                                              for i in range(ndims))});
    fpdtype_t tau;
    fpdtype_t nu_c = ${c['nu']};

% if spng_vis == 'mu':
    nu_c += spngmuin;
% endif

% for i in range(ndims):
    fout[${i}][${i + 1}] += -2.0*nu_c*(grad_uin[${i}][${i + 1}] - divu);
% endfor

% for i in range(ndims - 1):
    tau = -nu_c*(grad_uin[${i + 1}][1]
                                 + grad_uin[0][${i + 2}]);
    fout[${i + 1}][1] += tau;
    fout[0][${i + 2}] += tau;
% endfor

% if ndims == 3:
    tau = -nu_c*(grad_uin[1][3] + grad_uin[2][2]);
    fout[2][2] += tau;
    fout[1][3] += tau;
% endif
</%pyfr:macro>
