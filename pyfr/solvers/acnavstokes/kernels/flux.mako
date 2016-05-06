# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>
% for i in range(ndims):
    fout[${i}][${i+1}] += -2.0*${c['nu']}*(grad_uin[${i}][${i+1}]
                        - ${1.0/3.0}*(${' + '.join('grad_uin[{j}][{k}]'
                        .format(i=i,j=j,k=j+1) for j in range(ndims))}));
% endfor

% for i in range(ndims-1):
    fout[${i+1}][1] = fout[0][${i+2}] = -${c['nu']}*(grad_uin[${i+1}][1]
                                                    +grad_uin[0][${i+2}]);
% endfor

% if ndims == 3:
    fout[2][2] = fout[1][3] = -${c['nu']}*(grad_uin[1][3]+grad_uin[2][2]);
% endif
</%pyfr:macro>
