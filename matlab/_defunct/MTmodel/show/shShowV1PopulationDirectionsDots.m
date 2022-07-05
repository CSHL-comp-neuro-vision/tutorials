% shShowV1PopulationDirectionsDots(v1Neurons)
%
% Show the directions of the neurons in a V1 population as dots on a sphere
% in 3D Fourier space.
%
% Required arguments:
% v1Neurons         a matrix giving the parameters of V1 neurons in the
%                   standard format: each row gives the parameters for a
%                   different neuron; the first column gives the preferred
%                   direction of motion in radians with 0 = rightward; the
%                   second column gives the preferred ratio of temporal
%                   frequency to spatial frequency in cycles/frame and
%                   cycles/pixel, respectively.

function shShowV1PopulationDirectionsDots(v1Neurons)

v1Neurons(:,2) = atan3(v1Neurons(:,2), ones(size(v1Neurons,1), 1));
v1Neurons = sphere2rec(v1Neurons);
v1Neurons = v1Neurons*[1 0 0; 0 1 0; 0 0 -1];
v1Neurons = [v1Neurons; -v1Neurons];

[x,y,z] = sphere(50);
clf
axis([-1 1 -1 1 -1 1]);
hold on
surfl(x,y,z);
colormap(gray);
shading interp;
alpha(.6);
rotate3d on

for i = 1:size(v1Neurons, 1)
    plot3(v1Neurons(i, 2), v1Neurons(i, 1), v1Neurons(i, 3), 'ro', 'markerfacecolor', 'r');
end

axis off;
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
set(gca, 'units', 'pixels');
xlabel('x')
ylabel('y')
zlabel('t')
crossAxes('color', 'k');
labelCrossAxes('\omega_x', '\omega_y', '\omega_t', 'default', 'fontsize', 16);
rotate3d on
axis equal
axis vis3d
