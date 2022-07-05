function [pyr] = insertSpyrBand(pyr,pind,level,band,img)
% [LEV,IND] = insertSpyrBand(PYR,PYR_INDICES,LEVEL,BAND, IMG)
%
% Insert a band into a steerable pyramid.
%   LEVEL indicates the scale (1=finest).  BAND (optional) indicates 
%   which subband (1 = vertical, rest proceeding anti-clockwise).
%
% JED, 10/96

if (exist('band') ~= 1)
  band = 1;
end

nlevels=spyrHt(pind);
nbands=spyrNumBands(pind);

if ((level > nlevels) | (level < 1))
  error(sprintf('Bad level number: 	%d',level));
end

if ((band > nbands) | (band < 1))
  error(sprintf('Bad band number: 	%d',band));
end

firstband = 1 + band + nbands*(level-1);
pyr = insertPyrBand(pyr, pind, firstband,img);
return;

%%% Debug
ein=pgmRead('einstein.pgm');
[pyr,pind]=buildSpyr(ein,3,'sp3Filters');
showSpyr(pyr,pind);
pyr1=insertPyrBand(pyr,pind,6,zeros(size(pyrBand(pyr,pind,6))));
showSpyr(pyr1,pind);
pyr2=insertSpyrBand(pyr,pind,2,1,zeros(size(spyrBand(pyr,pind,2,1))));
showSpyr(pyr2,pind);
