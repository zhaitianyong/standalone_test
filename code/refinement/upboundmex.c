#include "mex.h"
#include <math.h>

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    double *X, *U, *err;
    int i, j, ij;
    double dist, temp, target;
    int nu, n , d;
    int ti, tj;
    
    X  = mxGetPr(prhs[0]);
    d  = mxGetM(prhs[0]);
    n  = mxGetN(prhs[0]);
    
    
    U = mxGetPr(prhs[1]);
    nu = mxGetM(prhs[1]);
    
    
       
    /* mxArray for output */
    plhs[0] = mxCreateDoubleMatrix(nu, 1, mxREAL);
    /* pointer to output data */
    err = mxGetPr(plhs[0]);
    
      
    /* lower bounds */
    for (i = 0 ; i < nu ; i++) {
        ti = (U[i] - 1) * d;
        tj = (U[i+nu] - 1) * d;
        target = U[i+2*nu];
        dist = 0;
        for (j = 0; j < d ; j++) {
            temp = X[ti+j] - X[tj+j];
            dist += temp * temp;
        }
        dist = sqrt(dist);
        temp = target - dist;
        if (temp < 0)
            err[i] = -temp;
    }
    
   
}

