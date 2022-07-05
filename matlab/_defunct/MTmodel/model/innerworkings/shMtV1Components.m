% v1Components = shMtV1Neuron(mtNeuron)     
%
% Find the parameters of V1 neurons that together feed forward to a given MT neuron.
%
% Required argument:
% mtNeuron          The parameters of an MT neuron. The first element is
%                   the preferred direction with 0 = right. The second is
%                   the preferred speed in pixels/frame.
%
% Output:
% v1Components      a 4x2 matrix. Each row contains the parameters of one
%                   of the V1 neurons that feeds forward to the MT neuron
%                   specified by mtNeuron. The first column gives the
%                   neuron's preferred direction in radians with 0 = right.
%                   The second column is the ratio of the neuron's
%                   preferred temporal frequency to spatial frequency in
%                   cycles/frame and cycles/pixel.

function v1Components = shMtV1Components(mtNeuron)

oldVel = mtNeuron;
mtNeuron(:, 1) = oldVel(:,2).*sin(oldVel(:,1));
mtNeuron(:, 2) = oldVel(:,2).*cos(oldVel(:,1));

if (sum(mtNeuron.^2) > 10*eps)
    d1 = [ -mtNeuron(1), -mtNeuron(2), sum(mtNeuron.^2)];
    d1 = d1/norm(d1);
    d2 = [ -mtNeuron(2), mtNeuron(1),  0          ];
    d2 = d2/norm(d2);
else
    d1 = [1 0 0];
    d2 = [0 1 0];
end
angles = pi*[0:3]'/(3+1);
v1Components = cos(angles)*d1 + sin(angles)*d2;

v1Components = rec2sphere(v1Components);
v1Components(:,3) = [];
v1Components(:,1) = mod(v1Components(:,1) + pi, 2*pi);
v1Components(:,2) = tan(v1Components(:,2));
