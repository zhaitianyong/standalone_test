#include "mex.h"
#include <math.h>

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    double *X, *E, *L, *U, *w, *f;
    int i, j, ij;
    double we, wl, wu, wr, *obj, dist, temp, target;
    int ne, nl, nu, n , d;
    double f_hb, f_tau, f_tal, f_vdw;
    
    int ti, tj, type;
    
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
    we = w[0];
    wl = w[1];
    wu = w[2];
    wr = w[3];
    
    f = mxGetPr(prhs[5]);
    
    f_hb  = f[0];
    f_tau = f[1];
    f_tal = f[2];
    f_vdw = f[3];
    
       
    /* mxArray for output */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    /* pointer to output data */
    obj = mxGetPr(plhs[0]);
    
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
        dist = sqrt(dist);
        temp = target - dist;
        obj[0] += we * temp * temp;
    }
    
    
    /* lower bounds */
    for (i = 0 ; i < nl ; i++) {
        ti = (L[i] - 1) * d;
        tj = (L[i+nl] - 1) * d;
        target = L[i+2*nl];
        type   = L[i+3*nl];
        dist = 0;
        for (j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        dist = sqrt(dist);
        temp = target - dist;
        if (temp > 0){
            if (type == 1){
                obj[0] += f_vdw * wl * temp * temp;
            }else if (type == -2){
                obj[0] += f_tal*wl * temp * temp;
            }
        }
    }
    
    
    /* upper bounds */
    for (i = 0 ; i < nu ; i++) {
        ti = (U[i] - 1) * d;
        tj = (U[i+nu] - 1) * d;
        target = U[i+2*nu];
        type = U[i+3*nu];
        dist = 0;
        for (j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        dist = sqrt(dist);
        temp = target - dist;
        if (temp < 0){
            if (type >= 0){
                obj[0] += wu * temp * temp;
            }else if (type == -1){
                obj[0] += f_hb*wu * temp * temp;
            }else if (type == -2){
                obj[0] += f_tau*wu * temp * temp;
            }
        }
    }
    
    
    /* regularization */
    for (i = 0 ; i < n ; i++) {
        ti = i*d;
        dist = 0;
        for (j = 0; j < d ; j++) {
            temp = X[ti+j];
            dist += temp * temp;
        }
        obj[0] += wr * dist;
    }
    
    
}

/*
 * %*****************************************************************
 * % Find the function value
 * % Originally by Dr. Kim Toh
 * %*****************************************************************
 * function objval = objfun(X,S,d_eq,D_lo,D_up,W_eq,W_lo,W_up,W_re)
 *
 * Xij     = X*S;
 * normXij = sqrt(sum(Xij.*Xij))';
 * objval  = W_eq*norm(normXij - d_eq)^2;
 *
 * %D = (distcalc(X)).^0.5;
 * %D = distcalc(X);
 * D = distmex(X);
 *
 * if ~isempty(D_lo)
 * [S_lo dd_lo] = lo_S_finder(D,D_lo);
 * if ~isempty(S_lo)
 * Xij = X*S_lo;
 * normXij = sqrt(sum(Xij.*Xij))';
 * objval = objval + W_lo*norm(normXij-dd_lo)^2;
 * end
 * end
 *
 * if isempty(D_up)
 * [S_up dd_up] = up_S_finder(D,D_up);
 * if ~isempty(S_up)
 * Xij = X*S_up;
 * normXij = sqrt(sum(Xij.*Xij))';
 * objval = objval + W_up*norm(normXij-dd_up)^2;
 * end
 * end
 *
 * objval = objval + W_re*sum(sum(X.^2));
 * %*****************************************************************
 */