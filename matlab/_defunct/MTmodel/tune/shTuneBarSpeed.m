% [xSpeed, yResponse] = shTuneBarSpeed(pars, neuron, stageName, nDataPoints,
%             minSpeed, maxSpeed, barDirection, barWidth, barEdgeWidth, rfRadius)


function [xSpeed, yResponse] = shtunebarspeed(varargin)

% The following arguments are optional and by default are 'default'
nDataPoints = 'default';
minSpeed = 'default';
maxSpeed = 'default';
barDirection = 'default';
barWidth = 'default';
barEdgeWidth = 'default';
rfRadius = 'default';

% parse the varargin
pars = varargin{1};
neuron = varargin{2};
stageName = varargin{3};
if nargin >= 4;     nDataPoints = varargin{4};          end;
if nargin >= 5;     minSpeed = varargin{5};             end;
if nargin >= 6;     maxSpeed = varargin{6};             end;
if nargin >= 7;     barDirection = varargin{7};         end;
if nargin >= 8;     barWidth = varargin{8};             end;
if nargin >= 9;     barEdgeWidth = varargin{9};         end;
if nargin >= 10;    rfRadius = varargin{10};              end;

% Assign default values where appropriate
if strcmp(nDataPoints, 'default');          nDataPoints = 7;                end;
if strcmp(minSpeed, 'default');             minSpeed = .25*neuron(2);       end;
if strcmp(maxSpeed, 'default');             maxSpeed = 4*neuron(2);         end;
if strcmp(barDirection, 'default');         barDirection = neuron(1);       end;
if strcmp(barWidth, 'default');             barWidth = 1;                   end;
if strcmp(barEdgeWidth, 'default');         barEdgeWidth = 11;               end; %4
if strcmp(rfRadius, 'default');             rfRadius = 20;                  end;


% END PARSING OF ARGUMENTS
targetStimSz = shGetDims(pars, stageName, [15 15 71]);
xSpeed = linspace(log2(minSpeed), log2(maxSpeed), nDataPoints);
xSpeed = 2.^xSpeed;
yResponse = zeros(1, length(xSpeed));
for i = 1:length(xSpeed)
    
    % generate the drifting bar stimulus
    nBarPasses = round(targetStimSz(3) .* xSpeed(i) ./ rfRadius);
    if nBarPasses == 0;
        nBarPasses = 1;
    end
    stimSz = [targetStimSz(1:2), ceil(nBarPasses*2*rfRadius/xSpeed(i))];
    thisBar = mkBar(stimSz, barDirection, xSpeed(i), barWidth, barEdgeWidth, 2*rfRadius, -rfRadius); %1*rfRadius
    thisBar = thisBar./20;
    
    % apply a lowpass filter to avoid aliasing.
    f = fftshift(fftn(ifftshift(thisBar)));
    fWindow = mkWedge(stimSz, [5 .4*stimSz(3)], neuron(1), pi/4, [5 .05*stimSz(3)]);
%     fWindow = mkWin(stimSz, [13 .65*stimSz(3)], [5 .05*stimSz(3)]);
    
%     
    f = f.*fWindow;
    thisBarFiltered = fftshift(ifftn(ifftshift(f)));
    thisBarFiltered = real(thisBarFiltered);

    % Pad the stimulus with zeros at the beginning and end.
    timePadSz = shGetDims(pars, 'mtPattern');
    timePadSz = timePadSz(3) - 1;
    thisStim = zeros(stimSz(1), stimSz(2), stimSz(3) + timePadSz);
    thisStim(:,:,timePadSz/2+1:end-timePadSz/2) = thisBar;
%     
%     saveFigure = gcf;    
%     figure(i);
%     clf reset
%     set(gcf, 'color', 'w');
%     f = fftshift(fftn(ifftshift(thisBar-mean2(thisBar))));
%     f = abs(f);
%     f = f./max2(f);
%     axis([-.5 .5 -.5 .5 -.5 .5]);
%     draw3dLevelSurfaces(f.*fWindow, 'lev', .1, 'alpha', .5', ...
%                         'color', [0 0 1], 'edgeColor', 'none');
%    
%     colorWts = shMtWts(neuron, pars);
%     hold on;
%     shShowV1NeuronSpectrum2(pars, pars.v1PopulationDirections, 'colorWts', ...
%                             colorWts, 'reso', 31, 'alpha', .5, 'lev', .8);
%     hold off;
%     axis square
%     axis off
%     axis vis3d
%     crossAxes('color', 'k');
%     labelCrossAxes('x', 'y', 't');
%     rotate3d on
%     drawnow;
%     
%     figure(saveFigure);
    
%     keyboard
    
    % Compute the neuron's response to this stimulus
    [pop, ind, res] = shModel(thisStim, pars, stageName, neuron);
    if strcmp(stageName, 'v1lin')
        res = sqrt(res.^2);
    end
%   yResponse(i) = mean2(shGetSubPop(res, ind, 1));
    yResponse(i) = mean2(shGetNeuron(res, ind, 1));

    
    
    % plot the responses so far
    semilogx(xSpeed(1:i), yResponse(1:i), 'r-', xSpeed(1:i), yResponse(1:i), 'k.');
    xlabel('speed (px/frame)'); ylabel('response');
    axis([min(xSpeed) max(xSpeed) min([0, yResponse]), max(.00001, 1.2*max(yResponse))]);
    drawnow
end