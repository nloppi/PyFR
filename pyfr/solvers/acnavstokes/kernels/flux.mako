# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
    fpdtype_t divu = ${1.0/3.0}*(${' + '.join('grad_uin[{i}][{j}]'
                                              .format(i=i,j=i + 1)
                                              for i in range(ndims))});
    fpdtype_t tau;

% for i in range(ndims):
    fout[${i}][${i + 1}] += ${-2.0*c['nu']}*(grad_uin[${i}][${i + 1}] - divu);
% endfor

% for i in range(ndims - 1):
    tau = ${-c['nu']}*(grad_uin[${i + 1}][1]
                                 + grad_uin[0][${i + 2}]);
    fout[${i + 1}][1] += tau;
    fout[0][${i + 2}] += tau;
% endfor

% if ndims == 3:
    tau = ${-c['nu']}*(grad_uin[1][3] + grad_uin[2][2]);
    fout[2][2] += tau;
    fout[1][3] += tau;
% endif
</%pyfr:macro>
