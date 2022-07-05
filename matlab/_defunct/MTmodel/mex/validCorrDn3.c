/* 
RES = validCorrDn3(IM, FILT, STEP);
  >>> See validCorrDn3.m for documentation <<<
  EPS, 1/98.
*/
#include <matrix.h>  /* Matlab matrices */
#include <mex.h>

#define notDblMtx(it) (!mxIsNumeric(it) || !mxIsDouble(it) || mxIsSparse(it) || mxIsComplex(it))

void mexFunction(int nlhs,	     /* Num return vals on lhs */
		 mxArray *plhs[],    /* Matrices on lhs      */
		 int nrhs,	     /* Num args on rhs    */
		 const mxArray *prhs[]     /* Matrices on rhs */
		 )
  {
  double *image,*filt, *result;
  int x_idim, y_idim, t_idim=1, n_idim=1;
  int x_fdim, y_fdim=1, t_fdim=1, n_fdim=1;
  int x_rdim, y_rdim, t_rdim, n_rdim;
  int x_step = 1;  /* default values */
  int y_step = 1;
  int t_step = 1;
  int rdims[4];
 
  mxArray *arg;
  double *mxMat;
  int dim;
  int const *pdim;
  if (nrhs<2) mexErrMsgTxt("requres at least 2 args.");

  /* ARG 1: IMAGE  */
  arg = prhs[0];
  if notDblMtx(arg) mexErrMsgTxt("IMAGE arg must be a non-sparse double float matrix.");
  image = mxGetPr(arg);
  
  dim =(int) mxGetNumberOfDimensions(arg);
  pdim = mxGetDimensions(arg);
  if ((dim < 2) || (dim > 4)) {
    mexErrMsgTxt("\n The IMAGE must have 2-4 dimensions! \n");
  }
 
  x_idim = (int) pdim[0]; /* X is inner index! */
  y_idim = (int) pdim[1];
  if (dim >= 3){
    t_idim = (int) pdim[2];
    if (dim == 4)
      n_idim = (int) pdim[3];
    }

/*  printf("%d, %d, %d, %d\n",x_idim,y_idim,t_idim,n_idim);*/

  /* ARG 2: FILTER */
  arg = prhs[1];
  if notDblMtx(arg) mexErrMsgTxt("FILTER arg must be non-sparse double float matrix.");
  filt = mxGetPr(arg);
 
  dim =(int) mxGetNumberOfDimensions(arg);
  pdim = mxGetDimensions(arg);
  x_fdim = (int) pdim[0]; 
  
  if (dim==2){
    y_fdim = (int) pdim[1];   
  } 
  if((dim==3) || (dim==4)){
    y_fdim = (int) pdim[1];  
    t_fdim = (int) pdim[2];
  }
  n_fdim = n_idim;
  
  if ((x_fdim > x_idim) || (y_fdim > y_idim) || (t_fdim > t_idim))
    {
    mexPrintf("Filter: [%d %d %d %d], Image: [%d %d %d %d]\n",
	      x_fdim,y_fdim,t_fdim, n_fdim,x_idim,y_idim,t_idim,n_idim);
    mexErrMsgTxt("FILTER dimensions larger than IMAGE dimensions.");
    }

  /* ARG 3 (optional): STEP */
  if (nrhs>2)
      {
      arg = prhs[2];
      if notDblMtx(arg) mexErrMsgTxt("STEP arg must be a double float matrix.");
      if (mxGetM(arg) * mxGetN(arg) != 3)
    	 mexErrMsgTxt("STEP vector must contain three elements.");
      mxMat = mxGetPr(arg);
      x_step = (int) mxMat[0];
      y_step = (int) mxMat[1];
      t_step = (int) mxMat[2];
      if ((x_step<1) || (y_step<1) || (t_step<1))
         mexErrMsgTxt("STEP values must be greater than zero.");
      }
	  
  x_rdim = (x_idim-x_fdim+x_step) / x_step;
  y_rdim = (y_idim-y_fdim+y_step) / y_step;
  t_rdim = (t_idim-t_fdim+t_step) / t_step;
  n_rdim = n_idim;
  rdims[0]=x_rdim; 
  rdims[1]=y_rdim;
  rdims[2]=t_rdim;
  rdims[3]=n_rdim;

  /*  mxFreeMatrix(plhs[0]); */

  plhs[0] = (mxArray *) mxCreateDoubleMatrix(x_rdim,y_rdim*t_rdim*n_rdim,mxREAL);
  mxSetDimensions(plhs[0],rdims,4);

  if (plhs[0] == NULL) mexErrMsgTxt("Cannot allocate result matrix");
  result = mxGetPr(plhs[0]);
 
  valid_filter(image, x_idim, y_idim, t_idim, n_idim, filt, x_fdim, 
	       y_fdim, t_fdim, x_step, y_step, t_step, result);

  return;
  } 



