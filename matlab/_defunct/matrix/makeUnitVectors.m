function uVectors = makeUnitVectors(theta)
% makeUnitVectors: returns 2d unit vectors, given an angle
%
% uVectors = makeUnitVectors(theta)

uVectors = [cos(theta); sin(theta)]';

