function Mout = translateAffine2(Min,x0,y0)
% Mout = translateAffine2(Min,x0,y0)
%
% Given an alignment matrix (Min, in homogeneous coords) computed
% from a section of an image that was offset from the upper
% left by (x0,y0), compute the alignment matrix (Mout) that should
% be applied to the original image.
%
% 8/13/97  dhb, eah  Wrote it.
% 8/31/97  dhb       Change name to translateAffine2.m

% Make offset vector in homogeneous coords
offsetVec0 = [x0-1 y0-1]';

% Compute offset
offsetVec = Min(1:2,1:2)*offsetVec0 - offsetVec0;

% Compute Mout
Mout = Min;
Mout(1,3) = Mout(1,3) - offsetVec(1);
Mout(2,3) = Mout(2,3) - offsetVec(2);
