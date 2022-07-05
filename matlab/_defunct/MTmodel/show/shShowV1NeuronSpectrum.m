% [f, preferredSin] = shShowV1NeuronSpectrum(pars, v1Neurons, levelToDraw,
%                           levelMode, fSz, multipleFilterStyle, patchArgs)
%
% Show the fourier spectrum of given model V1 neuron(s).
%
% Required arguments:
% pars          a standars shModel PARS structure
% v1Neurons     a matrix giving the paramters of the V1 neuron(s) whose
%               fourier spectrum you want to see, in the standard format:
%               Each row gives the parameters for a different neuron; the
%               first column gives the preferred direction in radians with
%               0 = right; the second row gives the preferred ratio of
%               temporal frequency to spatial frequency in cycles per frame
%               and cycles per pixel.
%
% Optional arguments:
% levelToDraw   The value of the function at which you want a level surface
%               drawn. DEFAULT = .5
% levelMode     if levelMode is 'norm', the levelToDraw level will be
%               relative to the maximum value of the function. Otherwise
%               the value will just be the absolute value of the function.
%               DEFAULT = 'norm'
% fSz           the size of the matrix containing the Fourier spectrum. The
%               larger fSz is, the higher the resolution of the final
%               image, but the longer the code will take to run.
% multipleFilterStyle       when you give the parameters of more than one
%               neuron in v1Neurons, their spectra are shown all at once.
%               If multipleFilterStyle is set to 'max', then at each point
%               in Fourier space, the maximum value of all the individual
%               spectra is used. If it is set to 'sum', the sum of the
%               individual spectra is used. 'max' results in a lumpy plot
%               that shows you where in fourier space each individual
%               spectrum lies; 'max' results in a smooth plot that shows
%               you the fourier spectrum of the population of neurons
%               together. DEFAULT = 'sum'
% patchArgs     any number of arguments that are passed onto the call to
%               PATCH that actually draws the spectra. 
%               DEFAULT = {'edgecolor', 'none'}
%
% Output:
% f             a matrix containing the values of the 3D fourier transform
%               of the filter(s);
% preferredSin  the parameters of the preferred sinusoid of the neuron(s)
%               given in v1Neurons. The first element is the preferred
%               direction of motion; the second element is the preferred
%               spatial frequency in cycles per pixel; the third element
%               is the preferred temporal frequency in cycles per frame.

function [f, preferredSin] = shShowV1NeuronSpectrum(varargin)

%%%%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%%%%
levelToDraw = 'default';
levelMode = 'default';
fSz = 'default';
multipleFilterStyle = 'default';
patchArgs = 'default';

pars = varargin{1};
v1Neurons = varargin{2};
if nargin >= 3;         levelToDraw = varargin{3};              end;
if nargin >= 4;         levelMode = varargin{4};                end;
if nargin >= 5;         fSz = varargin{5};               end;
if nargin >= 6;         multipleFilterStyle = varargin{6};      end;
if nargin >= 7;         patchArgs = varargin(7:end);            end;

if strcmp(levelToDraw, 'default');      levelToDraw = .5;           end;
if strcmp(levelMode, 'default');        levelMode = 'norm';         end;
if strcmp(fSz, 'default');       fSz = 51;            end;
if strcmp(multipleFilterStyle, 'default');  multipleFilterStyle = 'sum';        end;
if strcmp(patchArgs, 'default');        patchArgs = {};         end;

v1Neurons(:,1) = mod(v1Neurons(:,1)+pi, 2*pi);




%%%%%%%%%%%%%

tfilts = pars.v1TemporalFilters;
sfilts = pars.v1SpatialFilters;
fsz = size(tfilts, 1);

fs = zeros(10, fsz.^3);

pt = 0;
for ot = 0:3
    for ox = 0:3-ot
        oy = 3 - ox - ot;
        x = reshape(sfilts(:,ox+1), [1 fsz 1]);
        y = reshape(sfilts(:,oy+1), [fsz 1 1]);
        t = reshape(tfilts(:,ot+1), [1 1 fsz]);

        y = repmat(y, [1 fsz fsz]);
        x = repmat(x, [fsz 1 fsz]);
        t = repmat(t, [fsz fsz 1]);
        f = x.*y.*t;
        pt = pt+1;
        fs(pt, :) = f(:)';
    end
end


ws = shSwts(v1Neurons);
for i = 1:size(ws, 1)
    w = ws(i,:);
    f = w*fs;
    f = reshape(f, [fsz, fsz, fsz]);

    f = fftn(f, [fSz, fSz, fSz]);
    f = fftshift(f);
    f = f.*conj(f);
    if ~exist('Fs')
        Fs = zeros(size(f));
    end
    if strcmp(multipleFilterStyle, 'add')
        Fs = Fs + f;
    else
        Fs = max(Fs, f);
    end
end

if strcmp(levelMode, 'norm')
    f = Fs./max2(Fs);
else
    f = Fs;
end


clf
tmp = linspace(-.5, .5, fSz);
[x,y,z] = meshgrid(tmp, tmp, tmp);
if exist('levelToDraw') ~= 1
    levelToDraw = .5.*f./max2(f);
end
p = feval(@patch, isosurface(x,y,z,f,levelToDraw), 'edgecolor', 'none', patchArgs{:});
isonormals(x,y,z,f,p)
daspect([1 1 1]);
view(3);

axis off;
axis equal;
axis([-.5, .5, -.5, .5, -.5, .5]);
crossAxes('color', 'k', 'linewidth', 1);
labelCrossAxes('\omega_x', '\omega_y', '\omega_t', 'default', 'fontsize', 16);
camlight
lighting gouraud
set(gca, 'color', 'w');
set(gcf, 'color', 'w');
view(45, 45);
rotate3d on
axis vis3d


% calculate the point of maximum response
i= find(f == max2(f));
i = max2(i);
[i,j,k] = ind2sub([fSz,fSz,fSz], i);

xf = ((i - 1)./(fSz-1)).*1 - .5;
yf = ((j - 1)./(fSz-1)).*1 - .5;
tf = ((k - 1)./(fSz-1)).*1 - .5;
maxpos = [xf, yf, tf];


% generate the 3d sinusoid that best stimulates this neuron
S = zeros(fSz, fSz, fSz);
S(i, j, k) = fSz.^3;
cp = floor(fSz/2) + 1;
S(2*cp - i, 2*cp - j, 2*cp - k) = fSz.^3;
S = ifftshift(S);
preferredSin = ifftn(S);
preferredSin = real(preferredSin);
preferredSin = preferredSin - min2(preferredSin);
preferredSin = preferredSin./max2(preferredSin);
