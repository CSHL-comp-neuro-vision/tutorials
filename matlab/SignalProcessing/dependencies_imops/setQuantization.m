function res = setQuantization(in,nLevels)
% res = setQuantization(in,nLevels)
%
%AUTHOR:  Wandell
%DATE:    05.27.95
%PURPOSE:
%  Quantize the input function to a set of nLevels.
%  On return, the levels always run from 1 to nLevels.  You must
%  rescale them to get what you want.
%
%  This routine is like Rick's quantize (in imageCalc).  But it serves
%  a slightly different (and easier) purpose.
%

res = round(scale(in,1,nLevels));
