/* Objective Function Calculator
 * SPROS Refinment
 *
 * Babak Alipanahi
 * March 5, 2011
 * University of Waterloo
 */



#include "mex.h"
#include <math.h>

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    double *X, *E, *L, *U, *w, *G, *f;
    int i, j, ti, tj, type;
    double we, wl, wu, wr, dist, sdist, temp, target, tempG;
    double f_hb, f_tau, f_tal, f_vdw;

    int ne, nl, nu, n , d;
    
    X  = mxGetPr(prhs[0]);
    d  = mxGetM(prhs[0]);
    n  = mxGetN(prhs[0]);
    
    
    E = mxGetPr(prhs[1]);
    ne = mxGetM(prhs[1]);
    
    
    L = mxGetPr(prhs[2]);
    nl = mxGetM(prhs[2]);
    
    U = mxGetPr(prhs[3]);
    nu = mxGetM(prhs[3]);
    
    w  = mxGetPr(prhs[4]);
    we = 2*w[0];
    wl = 2*w[1];
    wu = 2*w[2];
    wr = 2*w[3];
    
    f = mxGetPr(prhs[5]);
    
    f_hb  = f[0];
    f_tau = f[1];
    f_tal = f[2];
    f_vdw = f[3];
    
    /* mxArray for output */
    plhs[0] = mxCreateDoubleMatrix(d, n, mxREAL);
    /* pointer to output data */
    G = mxGetPr(plhs[0]);
    
    
    /* equality constraints */
    for (i = 0 ; i < ne ; i++) {
        ti = (E[i] - 1) * d;
        tj = (E[i+ne] - 1) * d;
        target = E[i+2*ne];
        dist = 0;
        for (j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        sdist = sqrt(dist);
        tempG = (1 - target/sdist);
        for (j = 0; j < d ; j++) {
            G[ti+j] += we*(X[ti+j] - X[tj+j])*tempG;
            G[tj+j] -= we*(X[ti+j] - X[tj+j])*tempG;
        }
    }
    
    /* lower bounds */
    for (i = 0 ; i < nl ; i++) {
        ti = (L[i] - 1) * d;
        tj = (L[i+nl] - 1) * d;
        target = L[i+2*nl];
        type   = L[i+3*nl];
        dist   = 0;
        for (j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        sdist = sqrt(dist);
        tempG = (1 - target/sdist);
        if (sdist < target) {
            for (j = 0; j < d ; j++) {
                if (type == 1){
                    G[ti+j] += f_vdw*wl*(X[ti+j] - X[tj+j])*tempG;
                    G[tj+j] -= f_vdw*wl*(X[ti+j] - X[tj+j])*tempG;
                } else if (type == -2){
                    G[ti+j] += f_tal*wl*(X[ti+j] - X[tj+j])*tempG;
                    G[tj+j] -= f_tal*wl*(X[ti+j] - X[tj+j])*tempG;
                }else{
                    mexPrintf("Error! Unkonw Lower bound!\n");
                }
            }
        }
    }
    
    
    /* upper bounds */
    for (i = 0 ; i < nu ; i++) {
        ti = (U[i] - 1) * d;
        tj = (U[i+nu] - 1) * d;
        target = U[i+2*nu];
        type   = U[i+3*nu];
        dist   = 0;
        for (j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        sdist = sqrt(dist);
        tempG = (1 - target/sdist);
        if (sdist > target) {
            for (j = 0; j < d ; j++) {
                if (type >= 0){
                    G[ti+j] += wu*(X[ti+j] - X[tj+j])*tempG;
                    G[tj+j] -= wu*(X[ti+j] - X[tj+j])*tempG;
                } else if (type == -1){
                    G[ti+j] += f_hb*wu*(X[ti+j] - X[tj+j])*tempG;
                    G[tj+j] -= f_hb*wu*(X[ti+j] - X[tj+j])*tempG;
                } else if (type == -2){
                    G[ti+j] += f_tau*wu*(X[ti+j] - X[tj+j])*tempG;
                    G[tj+j] -= f_tau*wu*(X[ti+j] - X[tj+j])*tempG;
                }else{
                    mexPrintf("Error! Unkonw Upper bound!\n");
                }
            }
        }
    }
    
    
    /* regularization */
    for (i = 0 ; i < n ; i++) {
        ti = i*d;
        for (j = 0; j < d ; j++) {
            G[ti+j] += wr*X[ti+j];
        }
    }
    
    
    
    
}

