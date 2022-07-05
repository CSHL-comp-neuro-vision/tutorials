/* 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  File: svalidconvolve.c
;;;  Author: Eero Simoncelli
;;;  Description: 3D convolution code
;;;  Creation Date: Jan, 1998.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
*/

#include <stdio.h>
#include <math.h>
#include "convolve.h"

/* --------------------------------------------------------------------
  Correlate FILT with IMAGE, subsampling according to STEP parameters, 
  with values placed into RESULT array.  RESULT
  dimensions should be ceil((stop-start)/step).  

**** NEEDS TO BE MADE MORE EFFICIENT! ****
------------------------------------------------------------------------ */

int valid_filter(image, x_idim, y_idim, t_idim, n_idim,filt, x_fdim, 
		 y_fdim, t_fdim,x_step, y_step, t_step, result)
  register double *image;
  register int x_fdim, x_idim;
  register double *result;
  register int x_step, y_step, t_step;
  double *filt;
  int y_idim, y_fdim, t_fdim, t_idim;
  {
  register int x,y,t, base_image, pos_image;
  register int y_fil, x_fil,t_fil;
  int n_rdim = n_idim,n;
  register int x_rdim = (x_idim-x_fdim+x_step)/x_step;
  register int y_rdim = (y_idim-y_fdim+y_step)/y_step;
  register int t_rdim = (t_idim-t_fdim+t_step)/t_step;
  double sum;

  for(n=0; n< n_rdim; ++n){
    for(t = 0; t< t_rdim; ++t){
      for(y = 0; y < y_rdim; ++y){
	for(x = 0; x < x_rdim; ++x){
	  base_image = n*t_idim*y_idim*x_idim + 
       	               t*t_step*y_idim*x_idim + 
		       y*y_step*x_idim + 
		       x*x_step;
	  sum = 0.0;
	  for(t_fil=0; t_fil < t_fdim; ++t_fil){
	    for(y_fil=0; y_fil< y_fdim; ++y_fil){
	      pos_image = base_image + y_fil*x_idim + t_fil*x_idim*y_idim;
	      for(x_fil = 0; x_fil < x_fdim; ++x_fil,++pos_image){
		sum += image[pos_image] * filt[n*t_fdim*y_fdim*x_fdim +
					       t_fil*y_fdim*x_fdim +
					       y_fil*x_fdim +
					       x_fil];
		}	
	      }
	    }
	  result[n*t_rdim*y_rdim*x_rdim + t*y_rdim*x_rdim + y*x_rdim + x] = sum;
	  }
	}
      }
    }
  }
