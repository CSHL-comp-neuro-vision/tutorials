/* 
[] = Dsqr(MTX)
  Destructively replace all values in MTX by their squares.
  USE WITH EXTREME CAUTION: ARRAYS ARE OFTEN SHARED!!!

  EPS, 5/98.
*/

/* Matlab V4 types should be changed as follows:
  Matrix -> mxArray
  REAL -> mxREAL
  mxCreateFull ->  mxCreateDoubleMatrix
  */


#include <matrix.h>  /* Matlab matrices */
#include <mex.h>
#include <stddef.h>  /* NULL */
#define notDblMtx(it) (!mxIsNumeric(it) || !mxIsDouble(it) || mxIsSparse(it) || mxIsComplex(it))

void mexFunction(int nlhs,	     /* Num return vals on lhs */
		 mxArray *plhs[],    /* Matrices on lhs      */
		 int nrhs,	     /* Num args on rhs    */
		 const mxArray *prhs[]     /* Matrices on rhs */
		 )
  {
  double *mtx;
  int i, size;
  mxArray *arg;

  if (nrhs != 1) mexErrMsgTxt("requires 1 argument.");

  /* ARG 1: MATRIX  */
  arg = prhs[0];
  if notDblMtx(arg) mexErrMsgTxt("MTX arg must be a real non-sparse matrix.");
  mtx = mxGetPr(arg);
  size = (int) mxGetM(arg) * mxGetN(arg);

  for (i=0; i<size; i++)
      mtx[i] = mtx[i] * mtx[i];

  return;
  }      
