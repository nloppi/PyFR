# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

% if ndims == 2:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>

    // Velocity derivatives (grad[u,v])
    fpdtype_t u_x = grad_uin[0][1];
    fpdtype_t u_y = grad_uin[1][1];
    fpdtype_t v_x = grad_uin[0][2];
    fpdtype_t v_y = grad_uin[1][2];

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*${c['nu']}*(u_x - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_yy = -2*${c['nu']}*(v_y - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t t_xy = -${c['nu']}*(v_x + u_y);

    fout[0][1] += t_xx;     fout[1][1] += t_xy;
    fout[0][2] += t_xy;     fout[1][2] += t_yy;

</%pyfr:macro>

% elif ndims == 3:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout'>

    // Velocity derivatives (grad[u,v,w])
    fpdtype_t u_x = grad_uin[0][1];
    fpdtype_t u_y = grad_uin[1][1];
    fpdtype_t u_z = grad_uin[2][1];
    fpdtype_t v_x = grad_uin[0][2];
    fpdtype_t v_y = grad_uin[1][2];
    fpdtype_t v_z = grad_uin[2][2];
    fpdtype_t w_x = grad_uin[0][3];
    fpdtype_t w_y = grad_uin[1][3];
    fpdtype_t w_z = grad_uin[2][3];

    fpdtype_t nu_c = ${c['nu']};

    // Negated stress tensor elements
    fpdtype_t t_xx = -2*${c['nu']}*(u_x - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_yy = -2*${c['nu']}*(v_y - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_zz = -2*${c['nu']}*(w_z - ${1.0/3.0}*(u_x + v_y + w_z));
    fpdtype_t t_xy = -nu_c*(v_x + u_y);
    fpdtype_t t_xz = -nu_c*(u_z + w_x);
    fpdtype_t t_yz = -nu_c*(w_y + v_z);

    fout[0][1] += t_xx;     fout[1][1] += t_xy;     fout[2][1] += t_xz;
    fout[0][2] += t_xy;     fout[1][2] += t_yy;     fout[2][2] += t_yz;
    fout[0][3] += t_xz;     fout[1][3] += t_yz;     fout[2][3] += t_zz;

</%pyfr:macro>
% endif
