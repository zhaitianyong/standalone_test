#include "mex.h"
/*#include <math.h>*/
/* Distances
 * Computes the distances between each pair of n observations on a
 * k-dimensional
 * variable, x, using a quadratic distance metric defined by
 * dij = (xi-xj)(xi-xj)'
 * Paul Fackler
 * http://www.mathworks.com/matlabcentral/newsreader/view_thread/13163
 * some modifications by Babak Alipanahi
 * March 4, 2011
 */

void mexFunction(
        int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *x, *D, *xi, *xj, temp, *xend;
    int i, j, k, n, p, ij, np;
    
    if (!mxIsDouble(prhs[0]))
        mexErrMsgTxt("Input of inproper type");
    
    p = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    x = mxGetPr(prhs[0]);
    
    /* mxArray for output */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    /* pointer to output data */
    D = mxGetPr(plhs[0]); 
    
    xend = x + n*p;
    for (i = 0 ; i < n ; i++)
    {
        for (j = i; j < n ; j++)
        {
            ij = i + j*n;
            xi = x + i*p;
            xj = x + j*p;
            for (k = 0; k < p ; k++)
            {
                temp = *(xi + k) - *(xj + k);
                D[ij] += temp*temp;
            }
            D[j+n*i] = D[ij];
        }
    }
}
