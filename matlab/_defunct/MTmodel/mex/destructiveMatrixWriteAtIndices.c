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
  double *newValues;
  double *doubleStartIndex;
  int i, startIndex, size; 
  mxArray *arg;

  if (nrhs != 3) mexErrMsgTxt("requires 3 arguments.");

  /* ARG 1: MATRIX  */
  arg = prhs[0];
  if notDblMtx(arg) mexErrMsgTxt("MTX arg must be a real non-sparse matrix.");
  mtx = mxGetPr(arg);
  
  arg = prhs[1];
  if notDblMtx(arg) mexErrMsgTxt("MTX arg must be a real non-sparse matrix.");
  newValues = mxGetPr(arg);
  size = (int) mxGetM(arg) * mxGetN(arg);
  
  arg = prhs[2];
  if notDblMtx(arg) mexErrMsgTxt("MTX arg must be a real non-sparse matrix.");
  doubleStartIndex = mxGetPr(arg);
  startIndex = (int) doubleStartIndex[0];
  
  for (i=0; i<size; i++)
  	{  
      mtx[i+startIndex] = newValues[i];
  	}
  return;
  }