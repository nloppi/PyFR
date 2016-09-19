# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

void
axnpby(long *nrowp, long *ncolbp, long *ldimp,
       ${', '.join('fpdtype_t **xp' + str(i) for i in range(nv))},
       ${', '.join('double *ap' + str(i) for i in range(nv))})
{
    #define X_IDX_AOSOA(v, nv) ((ci/SOA_SZ*(nv) + (v))*SOA_SZ + cj)
    int ldim = *ldimp;

% for i in range(nv):
    fpdtype_t *x${i} = *xp${i};
    fpdtype_t  a${i} = *ap${i};
% endfor

    #pragma omp parallel
    {
        int align = PYFR_ALIGN_BYTES / sizeof(fpdtype_t);
        int rb, re, cb, ce, idx;
        fpdtype_t axn;
        loop_sched_2d(*nrowp, *ncolbp, align, &rb, &re, &cb, &ce);
        int nci = ((ce - cb) / SOA_SZ)*SOA_SZ;

        for (int r = rb; r < re; r++)
        {
            for (int ci = cb; ci < cb + nci; ci += SOA_SZ)
            {
                #pragma omp simd
                for (int cj = 0; cj < SOA_SZ; cj++)
                {
                % for k in subdims:
                    idx = r*ldim + X_IDX_AOSOA(${k}, ${ncola});
                    axn = ${pyfr.dot('a{l}', 'x{l}[idx]', l=(1, nv))};

                    if (a0 == 0.0)
                        x0[idx] = axn;
                    else if (a0 == 1.0)
                        x0[idx] += axn;
                    else
                        x0[idx] = a0*x0[idx] + axn;
                % endfor
                }
            }

            for (int ci = cb + nci, cj = 0; cj < ce - ci; cj++)
            {
            % for k in subdims:
                idx = r*ldim + X_IDX_AOSOA(${k}, ${ncola});
                axn = ${pyfr.dot('a{l}', 'x{l}[idx]', l=(1, nv))};

                if (a0 == 0.0)
                    x0[idx] = axn;
                else if (a0 == 1.0)
                    x0[idx] += axn;
                else
                    x0[idx] = a0*x0[idx] + axn;
            % endfor
            }
        }
    }
    #undef X_IDX_AOSOA
}
