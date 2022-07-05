% [v1Directions, dirsOverTime] = shParsdir(nNeurons, nSteps)       
%
% Get parameters for a population of V1 neurons that is evenly spread over
% sphere in 3D fourier space.
%
% There is no general solution for spreading n points over a sphere, so we
% use a dynamical model to place the points. The points all start at random
% points and then exert electrostatic and spring forces on each other. They
% are also damped by friction. Eventually they settle down to a final
% position.
%
% Required arguments:
% nNeurons          the number of neurons in the population
%
% Optional arguments:
% nSteps            the number of steps. The larger nSteps is, the longer
%                   program takes to run; but if nSteps is too small, the
%                   directions won't come close to a final position.
%                   DEFAULT = 400
% kPlane            the force exerted by the X-Y plane on the points on 
%                   the sphere. The larger kPlane is, the farther from the
%                   plane the points will lie. DEFAULT = .025
% kPole             the force extered by the Z-axis on the points.
%                   DEFAULT = 1
%
% Output:
% v1Directions      a matrix giving the parameters of the neurons in the V1
%                   population in the standard format: each row of
%                   v1Directions gives the parameters of a different neuron;
%                   the first column gives the preferred direction in
%                   radians with 0 = right; the second column gives the
%                   preferred ratio of temporal frequency to spatial
%                   frequency in cycles/frame and cycles/pixel,
%                   respectively.
% dirsOverTime      a matrix giving the directions of all the neurons at
%                   each time step. You can watch these as a movie using
%                   shShowV1PopulationDirectionsDotsMovie

function [v1Directions, dirsOverTime] = shParsv1dirs(varargin)

nSteps = 'default';
kPlane = 'default';
kPole = 'default';

nNeurons = varargin{1};
if nargin >= 2;         nSteps = varargin{2};       end
if nargin >= 3;         kPlane = varargin{3};       end
if nargin >= 4;         kPole = varargin{4};        end

if strcmp(nSteps, 'default');       nSteps = 400;       end
if strcmp(kPlane, 'default');       kPlane = .025;      end
if strcmp(kPole, 'default');        kPole = 1;          end

% END PARSING OF ARGUMENTS


nNeurons = nNeurons.*2;

theta0 = rand(nNeurons, 1);
theta0 = theta0 - .5;
theta0 = theta0.*2*pi;

phi0 = rand(nNeurons, 1);
phi0 = phi0 - .5;
phi0 = phi0.*pi;

x0 = sphere2rec([theta0, phi0]);
x0 = x0';

v0 = zeros(3,nNeurons);

xnull = 2;
ks = 1;   % 1
kd = .6;  % .6
h = .1;   % .1

[dirsOverTime, vs, ts, s, c, d] = shParsdirsrk(x0, v0, xnull, ks, kd, kPole, kPlane, h, nSteps);


% the directions are returned in cartestian coordinates. Convert them into
% the special coordinate system used in the model.
v1Directions = squeeze(dirsOverTime(:,1:nNeurons/2,end))';
v1Directions = rec2sphere(v1Directions);
v1Directions(:,1) = mod(v1Directions(:,1) + pi, 2.*pi);
v1Directions = v1Directions(:, 1:2);

w = find(v1Directions(:,2) < 0);
v1Directions(w,1) = mod(v1Directions(w,1) + pi, 2*pi);
v1Directions(w,2) = -1.*v1Directions(w,2);
v1Directions(:,2) = tan(v1Directions(:,2));










%%%%%%%% HERE'S WHERE THE ACTUAL RUNGA-KUTTA (SP?) GOES ON!
function [dirsOverTime,vs,ts, springs, charges, drags, planes, pole] = shParsdirsrk(x, v, x0, k, kd, kPole, kPlane, h, steps)

N = size(x, 2);
dirsOverTime = zeros(3, N, steps);
vs = zeros(3, N, steps);
ts = zeros(steps);
dirsOverTime(:,:,1) = x;
vs(:,:,1) = v;
ts(1) = 0;

for nNeurons = 1:steps-1
    t = ts(nNeurons);
    x = dirsOverTime(:,:,nNeurons);
    v = vs(:,:,nNeurons);

    k1 = rkv(t, x, v);
    [j1, s1, c1, d1] = rka(t, x, v, x0, k, kd, kPole, kPlane);

    k2 = rkv(t+h/2, x+(h/2)*k1, v+(h/2)*j1);
    [j2, s2, c2, d2] = rka(t+h/2, x+(h/2)*k1, v+(h/2)*j1, x0, k, kd, kPole, kPlane);

    k3 = rkv(t+h/2, x+(h/2)*k2, v+(h/2)*j2);
    [j3, s3, c3, d3] = rka(t+h/2, x+(h/2)*k2, v+(h/2)*j2, x0, k, kd, kPole, kPlane);

    k4 = rkv(t+h, x+(h/2)*k2, v+(h/2)*j2);
    [j4, s4, c4, d4] = rka(t+h, x+(h/2)*k2, v+(h/2)*j2, x0, k, kd, kPole, kPlane);

    x = x + h*k1/6 + h*k2/3 + h*k3/3 + h*k4/6;
    v = v + h*j1/6 + h*j2/3 + h*j3/3 + h*j4/6;
    t = t + h;

    springs(:,:,nNeurons+1) = s1/6 + s2/3 + s3/3 + s4/6;
    charges(:,:,nNeurons+1) = c1/6 + c2/3 + c3/3 + c4/6;
    drags(:,:,nNeurons+1) = d1/6 + d2/3 + d3/3 + d4/6;

    rxy = sqrt(sum(x(1:2,:).^2));
    phi = atan(x(3,:)./rxy);
    theta = x(2,:)./x(1,:);
    theta = atan(theta);

    r = sqrt(sum(x.^2));
    r(r>1) = 1;

    newx = [r.*cos(theta).*cos(phi); r.*sin(theta).*cos(phi); r.*sin(phi)];
    w = (sign(x) ~= sign(newx));
    newx(w) = newx(w)*-1;
    x = newx;

    % match pairs; first N/2 correspond to second N/2
    xm = x(:,1:N/2) + x(:,N/2+1:end);
    x = x - repmat(xm./2, [1 2]);

    vm = v(:,1:N/2) + v(:,N/2+1:end);
    v = v - repmat(vm./2, [1 2]);

    dirsOverTime(:,:,nNeurons+1) = x;
    vs(:,:,nNeurons+1) = v;
    ts(nNeurons+1) = t;
end


function v = rkv(t,x,v)
v = v;


function [a, spring, charge, drag, plane, pole] = rka(t,x,v,x0,k,kd, kPole, kPlane)

% first, calculate the distances between the particles
N = size(x,2);
xd = zeros(N,N,3);
vd = zeros(N,N,3);
for i = 1:N
    for j = 1:N
        xd(i,j,:) = x(:,i) - x(:,j);
        vd(i,j,:) = v(:,i) - v(:,j);
    end
end

xdm = xd.^2;
xdm = sum(xdm,3);
xdm = sqrt(xdm);
xdm = repmat(xdm, [1 1 3]);

vdm = vd.^2;
vdm = sum(vdm, 3);
vdm = sqrt(vdm);
vdm = repmat(vdm, [1 1 3]);             % where is VDM used???

warning off MATLAB:divideByZero
xd = xd./xdm;
w = isnan(xd);
xd(w) = 0;



spring = xdm - repmat(x0, size(xdm));
% spring = spring + kd.*vdm;            % VDM used here.
spring = -k.*(spring);
spring = spring.*xd;
spring(w) = 0;
spring = sum(spring, 2);
spring = reshape(spring, [N,3])';

% charge = -k.*xdm.*xd;
% charge(w) = 0;
% charge = sum(charge, 2);
% charge = reshape(charge, [N,3])';

charge = zeros(size(spring));
%%%%%

vm = sqrt(sum(v.^2));
vm = repmat(vm, [3, 1]);
vnorm = v./vm;
drag = -kd.*vnorm;
drag(isnan(drag)) = 0;


%%%% pole %%%%


pole = x;
pole(3,:) = 0;
fpole = sqrt(sum(pole.^2));
pole = repmat((kPole./(fpole.^2)), [3 1]) .* pole./repmat(fpole, [3 1]);


%%%% PLANE %%%%

plane = x;
plane(1:2,:) = 0;
fplane = sqrt(sum(plane.^2));
plane = repmat((kPlane./(fplane.^2)), [3 1]) .* sign(plane);


a = spring + charge + drag + plane + pole;

