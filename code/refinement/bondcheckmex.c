#include "mex.h"
#include <math.h>

void absolute(double *in)
{
    if (in[0] < 0)
        in[0] = -in[0];
}


void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    double *X, *E, *err;
    int i, j, ij;
    double dist, temp, target;
    int ne, n , d;
    int ti, tj;
    
    X  = mxGetPr(prhs[0]);
    d  = mxGetM(prhs[0]);
    n  = mxGetN(prhs[0]);
    
    
    E = mxGetPr(prhs[1]);
    ne = mxGetM(prhs[1]);
    
    
       
    /* mxArray for output */
    plhs[0] = mxCreateDoubleMatrix(ne, 1, mxREAL);
    /* pointer to output data */
    err = mxGetPr(plhs[0]);
    
      
    /* lower bounds */
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
        absolute(&temp);
        if (temp > 0.01)
            err[i] = temp;
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