% theta = atan3(y, x)      Get the four quadrant atan with 0 <= theta <= 2.*pi
%
% atan3 is slightly different than atan2 in that the output ranges from 0
% to 2*pi instead of from -pi to pi.
%
% SEE ALSO: ATAN, ATAN2

function theta = atan3(y, x)

theta = atan2(y,x);
theta(theta<0) = theta(theta<0) + 2*pi;