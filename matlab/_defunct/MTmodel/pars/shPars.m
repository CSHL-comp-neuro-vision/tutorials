% pars = shPars         get default parameters structure for the SH model.
%
% Feel free to create a new function using this as a template.
%
% Choosing certain parameters outside of certain ranges will lead to
% crashes:
%
% pars.v1C50 should be below around .6, though depending on other
% parameters it may be safe to push it higher.
%
% pars.mtC50 should be equal to pars.v1C50. It is possible in some cases to
% move it away from pars.v1C50, but this usually leads to unusual model
% behavior.
%
% pars.mtBaseline must be between .001 and 1. Values close to 1 lead to
% very strange model behavior.
%
% SEE ALSO: shParsScaleFactors, shParsV1PopulationDirections

function pars = shPars

% load some of the paramters that are big matrices that are no fun to type
% into this file when you change them.
directoryContainingThisFile = which('shPars');
w = find(directoryContainingThisFile == '/');
directoryContainingThisFile = directoryContainingThisFile(1:w(end));
load([directoryContainingThisFile, 'defaultParameters.mat']);

%%%%% NOW WE GET STARTED: V1
pars.nScales = 1;
pars.v1SpatialFilters = v1SpatialFilters;       % Linear filters used to compute V1 responses. Stored in defaultParameters.mat
pars.v1TemporalFilters = v1TemporalFilters;     % Linear filters used to compute V1 responses. Stored in defaultParameters.mat
pars.v1PopulationDirections = v1PopulationDirections;       % Parameters for neurons in the V1 population. Stored in defaultParameters.mat
pars.v1Baseline = 0;                            % Additive constant in V1. Always 0. Included out of fidelity to the original paper.
pars.v1ComplexFilter = mkGaussianFilter(1.6);   % Blurring filter used to make complex cell responses phase invariant.
pars.v1NormalizationType = 'tuned';             % Choices: 'tuned', 'untuned', and 'off';
                                                % 'untuned' and 'off' are diagnostic settings and shouldn't be used unless you know what you're about.
pars.v1NormalizationSpatialFilter = mkGaussianFilter(-1);   % Blurring filter to make the normalization pool larger than the CRF.
pars.v1NormalizationTemporalFilter = mkGaussianFilter(-1);  % Blurring filter to make the normalization signal pool over time.
pars.v1C50 = .1;                                % Contrast at which V1 neurons have half maximal response to a drifting grating.

%%%%% AND ON TO MT
pars.mtPopulationVelocities = mtPopulationVelocities;   % Preferred velocities of neurons in the MT population; stored in defaultParameters.mat
pars.mtSpatialPoolingBeforeThreshold = 1;           % Is spatial pooling performed before or after the half wave rectification of MT responses?
pars.mtSpatialPoolingFilter = mkGaussianFilter(3);  % MT spatial pooling filter
pars.mtNormalizationType = 'tuned';             % choices are 'tuned' and 'self'. 'self' is currently a diagnostic setting for those who know what they're about.
pars.mtNormalizationSpatialFilter = mkGaussianFilter(-1);   % Filter for spatial pooling of the MT normalization signal.
pars.mtNormalizationTemporalFilter = mkGaussianFilter(-1);  % Filter for temporal pooling of the MT normalization signal.
pars.mtC50 = pars.v1C50;                        % Contrast at which MT neurons have half maximum response to full field drifting gratings.
                                                % Model is unstable if v1C50 ~= mtC50.
pars.mtBaseline = .1;                           % Baseline response of MT neurons.
pars.mtExponent = 2;                            % Exponent to which MT neuron responses are raised.

%%%% COMPUTE SCALE FACTORS
pars = shParsScaleFactors(pars);
