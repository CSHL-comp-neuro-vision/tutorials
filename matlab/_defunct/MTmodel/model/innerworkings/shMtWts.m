% wts = shMtWts(mtNeurons, pars)    
%
% Compute the weights on V1 population responses for given MT neurons.
%
% Required arguments:
% mtNeurons         the parameters of the MT neurons for which you want the
%                   weights, in the standard format.
% pars              a parameters structure.
%
% Output:
% wts               wts on the responses of neurons in the V1 population

function wts = shMtWts(mtNeurons, pars)

wts = [];
for i = 1:size(mtNeurons, 1)
    dirs = shMtV1Components(mtNeurons(i,:));
    tmp = sum(shQwts(dirs) * pinv(shQwts(pars.v1PopulationDirections)));
    tmp = tmp-mean(tmp);
    wts = [wts; tmp];
end
