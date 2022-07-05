function [R,S,E,A,G,I] = normalizationModel(p,stim,attend)

%[R,S,E,A,G,I] = normalizationModel(p,stim,[attend])

%(S) Stimulus 
S = makeNeuralImage(p,stim);

%(E) Convolve stimulus with RF to get 'stimulus drive' 
E = convolveImage(p,S,p.e.x.width,p.e.theta.width);

%(A) Attention field 
if exist('attend','var')
    A = makeNeuralImage(p,attend);
else %No attention
    A = ones(size(E));  
end

%(G) Multiply drive by attention field to get the 'excitatory drive' 
G = E.*A;

%(I) Convolve excitatory drive by inhibitory filters for 'suppressive drive' 
I = convolveImage(p,G,p.i.x.width,p.i.theta.width);

%(R) Pass through normalization stage to get the population response
R = G./(I+p.sigma);


