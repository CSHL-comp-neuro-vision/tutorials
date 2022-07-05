% Image Registration Tools.
% Last update:   18-Dec-1997
% Author: David Heeger, heeger@stanford.edu
% 
% See README file for brief description.
% Type "help <command-name>" to get help on individual commands.
%
% ---------------------------------------------------------------
% circularShift3 - shift a 3D array by integer number of samples.
% computeDerivatives2 - internal function that computes x-, y-,
%			   and t-derivatives of a pair of images.
% computeDerivatives3 - internal function that computes x-, y-,
%			   z-, and t-derivatives of a pair of
%			   volumes.
% convXYsep - 2D convolution on each sub-image of a volume.
% convZ - 1D convolution in the Z direction on a volume.
% 
% estTrans3 - estimate translation motion on a pair of volumes.
% estAffine2 - estimate 2D affine registration on a pair of images. 
% estAffine3 - estimate 3D affine registration on a pair of volumes.
% estAffineIter2 - iterative estimate of 2D affine registration
%			   on a pair of images.
% estAffineIter3 - iterative estimate of 2D affine registration
%			   on a pair of volumes.
% estAffineMulti2 - multiscale/iterative estimate of 2D affine
%			   registration on a pair of images.
% estAffineMulti3 - multiscale/iterative estimate of 2D affine
%			   registration on a pair of volumes. 
% 
% reduce - separable convolution and subsampling by a factor of 2. 
% translateAffine2 - Given an 2D affine transform estimated
%                          from a section of an image, compute
%                          the alignment matrix that should be
%                          applied to the original entire image. 
% warpAffine2 - warp an image according to a 2D affine transform.
% warpAffine3 - warp a volume according to a 3D affine transform.
