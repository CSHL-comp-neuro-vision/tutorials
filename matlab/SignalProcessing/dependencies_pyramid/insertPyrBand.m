function [pyr] = insertPyrBand(pyr,pind,band,img)
% insertPyrBand(pyr,pind,band,img)
%
% Insert an image into a subband of a pyramid
%    (gaussian, laplacian, wavelet, or steerable)
%
% JED, 10/96

ind = pyrBandIndices(pind,band);
target = reshape(img,pind(band,1) * pind(band,2),1);

pyr(ind) = target;
return;

%%% Debug
ein=pgmRead('einstein.pgm');
[pyr,pind]=buildLpyr(ein);
showLpyr(pyr,pind);
pyr1=insertPyrBand(pyr,pind,2,zeros(size(pyrBand(pyr,pind,2))));
showLpyr(pyr1,pind);
