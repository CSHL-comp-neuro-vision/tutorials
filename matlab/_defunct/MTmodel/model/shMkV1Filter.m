% [filter, S] = shMkV1Filter(pars, v1Neurons)
%
% Make the linear filter that is the front end of a given model V1 neuron.
%
% N.B.: The shModel code actually never calls shMkV1Filter, although it
% uses equivalent code. This function is provided simply to allow you to
% inspect the filters visually.
%
% Required arguments:
% pars          a standard shModel PARS structure.
% v1Neurons     a matrix containing the neurons whose linear filters you 
%               want to inspect, in standard form. Each row gives the
%               parameters for a different neuron; the first column gives
%               the preferred direction in radians with 0 = right; the
%               second column gives the preferred ratio of temporal
%               frequency to spatial frequency in cycles/frame and
%               cycles/pixel, respectively.
% 
% Outputs:
% filters       a 4D matrix containing the filters. Use
%               squeeze(filters(:,:,:,n)) to get the filter for the neuron
%               given in the nth row of v1Neurons. If v1Neurons contains
%               just one row, filters will just be a 3D matrix containing
%               the appropriate filter.
% S             A 4D matrix containing the separable filters that are
%               combined to make the V1 neuron filters.
%
% Example of use:
% pars = shPars;
% filters = shMkV1Filter(pars, [0 1; 0 2]);
% flipBook(squeeze(filters(:,:,:,1)), 'default', .1);
% flipBook(squeeze(filters(:,:,:,2)), 'default', .1);
        
function [filters, S] = shMkV1Filter(pars, v1Neurons)

sfilts = pars.v1SpatialFilters;
tfilts = pars.v1TemporalFilters;

order = 3;
fsz = size(sfilts, 1);
n = 1;
S = zeros(fsz, fsz, fsz, 10);


for torder = 0:order
    for xorder = 0:(order-torder)
        yorder = order - torder - xorder;
        tfilt = flipud(tfilts(:,torder+1));
        xfilt = sfilts(:,xorder+1);
        yfilt = flipud(sfilts(:,yorder+1));

        tmp = yfilt * xfilt';
        %         keyboard
        for i = 1:size(tfilt, 1)
            S(:,:,i,n) = tmp .* tfilt(i);      % store the result in S
        end
        n = n+1;
    end
end

S = reshape(S, [fsz.^3, 10])';
filters = shSwts(v1Neurons) * S;
filters = reshape(filters', [fsz, fsz, fsz, size(v1Neurons, 1)]);
S = reshape(S', [fsz, fsz, fsz, 10]);


flipBook(filters);