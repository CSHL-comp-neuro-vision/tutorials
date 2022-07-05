% pars = shParsSetScaleFactors(pars)
%
% Set scale factors for the pars structure and pick pars.mtalpha
%
% shParsSetScaleFactors chooses the pars.scaleFactors values so that every
% neuron responds to its preferred grating (full field at full contrast)
% with a response of 1. 
%
% It also chooses parameters so that pars.v1C50 and pars.mtC50 are actually
% the contrasts at which the response of neurons is at half maximum.

function pars = shParsSetScaleFactors(pars)

stagename = 'mtpattern';
dims = shGetDims(pars, stagename, [1 1 31]);
pref = mt2sin([0 1]);
s = mkSin(dims, pref(1), pref(2), pref(3), 1);
s = (s+1)./2;

pars.scaleFactors.v1Linear = 1;
pars.scaleFactors.v1FullWaveRectified = 1;
pars.scaleFactors.v1Blur = 1;
pars.scaleFactors.v1NormalizationPopulationK = 1;
pars.scaleFactors.v1NormalizationStrength = 1;
pars.scaleFactors.v1Complex = 1;
pars.scaleFactors.mtLinear = 1;
pars.scaleFactors.mtHalfWaveRectification = 1;
pars.scaleFactors.mtNormalization = 1;
pars.scaleFactors.mtNormalizationStrength = 1;
pars.mtAlpha = 0;

[pop, ind, res] = shModel(s, pars, 'v1lin', [0 1]);
pars.scaleFactors.v1Linear = 1./max(shGetNeuron(res, ind));

[pop, ind, res] = shModel(s, pars, 'v1fullrect', [0 1]);
pars.scaleFactors.v1FullWaveRectified = 1./mean(shGetNeuron(res, ind));

[pop, ind, res] = shModel(s, pars, 'v1blur', [0 1]);
pars.scaleFactors.v1Blur = 1./mean(shGetNeuron(res, ind));

[pop, ind, nume, deno, res, resnume, resdeno] = shModel(s, pars, 'v1complex', [0 1]);
pars.scaleFactors.v1NormalizationPopulationK = 1./mean(shGetNeuron(resdeno, ind));
pars.scaleFactors.v1NormalizationStrength = 1 - 2.*pars.v1C50.^2;
pars.scaleFactors.v1Complex = 1 - pars.v1C50.^2;

[pop, ind, res] = shModel(s, pars, 'mtlin', [0 1]);
pars.scaleFactors.mtLinear = 1./mean(shGetNeuron(res, ind));

[pop, ind, res] = shModel(s, pars, 'mtlin', [0 1]);
pars.temporaryPopVec = mean(shGetNeuron(pop, ind), 2);
x = linspace(.0001, 10, 10000);
y = zeros(size(x));
y(1) = shParsRhoError(x(1), pars);
y(2) = shParsRhoError(x(2), pars);
i = 3;
while i < length(x)
    y(i) = shParsRhoError(x(i), pars);
    if y(i-1) < y(i-2) & y(i-1) < y(i) & y(i-1) < .01
        pars.mtAlpha = x(i-1);
        i = length(x)+10;
    end
    i = i+1;
end
[junk, M, N] = shParsRhoError(pars.mtAlpha, pars);
pars = rmfield(pars, 'temporaryPopVec');

pars.scaleFactors.mtHalfWaveRectification = 1./((1+pars.mtAlpha).^pars.mtExponent);
pars.scaleFactors.mtPattern = M./pars.scaleFactors.mtHalfWaveRectification;
pars.scaleFactors.mtNormalizationStrength = N./pars.scaleFactors.mtHalfWaveRectification;



%%%%%%
% pars.temporaryPopVec = mean(shGetNeuron(pop, ind), 2);
% x = linspace(.0001, 20, 10000);
% y = zeros(size(x));
% for i = 1:length(x);
%     y(i) = shParsRhoError(x(i), pars);
% end
% keyboard



function [res, M, N] = shParsRhoError(rho, pars)

b = pars.mtBaseline;
mtExponent = pars.mtExponent;
mt50 = pars.mtC50;
v150 = pars.v1C50;
popVec = pars.temporaryPopVec;
nPopulation = length(pars.mtPopulationVelocities);

%%%% CALCULATE V1 RESPONSE

v1ResponseNumerator = (1-v150.^2).*mt50.^2;
v1ResponseDenominator = (1-2.*v150.^2).*mt50.^2 + v150.^2;

%%%% CALCULATE M & N

nastySumOne = popVec + rho;
nastySumOne(nastySumOne < 0) = 0;
nastySumOne = sum(nastySumOne.^mtExponent);

mNumerator = b.*rho.^(-mtExponent).*mt50.^2.*(-nastySumOne+nPopulation.*rho.^mtExponent);
mDenominator = nPopulation.*b.*(1+rho).^mtExponent - nastySumOne;
M = mNumerator./mDenominator;

nNumerator = rho.^(-mtExponent).*mt50.^2.*(rho.^mtExponent-b.*(1+rho).^mtExponent);
nDenominator = nPopulation.*b.*(1+rho).^mtExponent - nastySumOne;
N = nNumerator./nDenominator;

%%%% NOW GET THE FULL EQUATION

nastySum = popVec.*(v1ResponseNumerator./v1ResponseDenominator) + rho;
nastySum(nastySum < 0) = 0;
nastySum = sum(nastySum.^mtExponent);

mtResponseNumerator = M.*(v1ResponseNumerator./v1ResponseDenominator + rho).^mtExponent;
mtResponseDenominator = N.*nastySum + mt50.^2;

res = abs(((1-b)./2 + b) - mtResponseNumerator./mtResponseDenominator);
