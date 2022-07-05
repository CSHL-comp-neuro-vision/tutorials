
function updatePath(save)
% 
% function updatePath(save)
% 
% updates Matlab's path for CSHL tutorials.
%   INPUT
%       save    0, the updated path will not be saved (default). you need to run the
%                  function again the next time you start Matlab. 
%               1, save the updated path permanently. 
%               
% RK, 6/2010
% GMB,6/18/10 added default case of no input variables.
% 

if ~exist('save','var')
    save = 0;
end

d = fileparts(which('updatePath.m'));
addpath(genpath(d));
if save
    savepath;
end


