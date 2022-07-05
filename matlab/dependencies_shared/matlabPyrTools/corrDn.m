% RES = corrDn(IM, FILT, EDGES, STEP, START, STOP)
%
% Compute correlation of matrices IM with FILT, followed by
% downsampling.  These arguments should be 1D or 2D matrices, and IM
% must be larger (in both dimensions) than FILT.  The origin of filt
% is assumed to be floor(size(filt)/2)+1.
%
% EDGES is a string determining boundary handling:
%    'circular' - Circular convolution
%    'reflect1' - Reflect about the edge pixels
%    'reflect2' - Reflect, doubling the edge pixels
%    'repeat'   - Repeat the edge pixels
%    'zero'     - Assume values of zero outside image boundary
%    'extend'   - Reflect and invert
%    'dont-compute' - Zero output when filter overhangs input boundaries
%
% Downsampling factors are determined by STEP (optional, default=[1 1]),
% which should be a 2-vector [y,x].
%
% The window over which the convolution occurs is specfied by START
% (optional, default=[1,1], and STOP (optional, default=size(IM)).
%
% NOTE: this operation corresponds to multiplication of a signal
% vector by a matrix whose rows contain copies of the FILT shifted by
% multiples of STEP.  See upConv.m for the operation corresponding to
% the transpose of this matrix.

% Eero Simoncelli, 6/96, revised 2/97.

% Geoff Boynton 6/10, edited to use 'convolve2' instead of mex function 'rconv2'.
function res = corrDn(im, filt, edges, step, start, stop)

%% NOTE: THIS CODE IS NOT ACTUALLY USED! (MEX FILE IS CALLED INSTEAD)

%fprintf(1,'WARNING: You should compile the MEX version of "corrDn.c" found in the MEX subdirectory of matlabPyrTools, and put it in your matlab path.  It is MUCH faster, and provides more boundary-handling options.\n');

%------------------------------------------------------------
%% OPTIONAL ARGS:

if (exist('edges') == 1)
    if (strcmp(edges,'reflect1') ~= 1)
        % warning('Using REFLECT1 edge-handling (use MEX code for other options).');
    end
end

if (exist('step') ~= 1)
    step = [1,1];
end

if (exist('start') ~= 1)
    start = [1,1];
end

if (exist('stop') ~= 1)
    stop = size(im);
end

%------------------------------------------------------------

% Reverse order of taps in filt, to do correlation instead of convolution
filt = filt(size(filt,1):-1:1,size(filt,2):-1:1);
switch edges
    case 'circular'
        tmp = convolve2(im,filt,'wrap');
    case 'reflect1' | 'reflect2'
        tmp = convolve2(im,filt,'reflect');
    case 'dont_compute'
        tmp = convolve2(im,filt,'same');
    otherwise
        tmp = convolve2(im,filt,'same');
end

res = tmp(start(1):step(1):stop(1),start(2):step(2):stop(2));
