% dims = shGetDims(pars, stageName, outputDims)   
%
% Find the size your stimulus must be to successfully run the MT model.
%
% Because the successive stages of the model perform repeated edge-dropping
% convolutions on the stimulus, your stimulus will be bigger than the
% output of the model. If your stimulus is too small, the program will crash
% when it runs out of dimensions to trim away. This program figures out
% how big your stimulus will have to be to get the job done.
%
% Required arguments:
% pars          a standard MT model parameters file.
% stageName     the name of the stage whose output you wish to compute.
%               Standard choices: 'v1Complex', 'mtPattern'.
%
% Optional arguments:
% outputDims    the size of the output you want from shModel in [Y X T]
%               coordinates. DEFAULT = [1 1 1].
%
% Output:
% dims          a three vector specifying the required stimulus size in 
%               [Y X T] coordinates.



function dims = shGetDims(pars, stageName, outputDims)

if exist('stage') ~= 1
    stage = 'mtnorm';
end

if exist('outputDims') ~= 1
    outputDims = [1 1 1];
end

dims = outputDims;

switch lower(stageName)
    case 'v1lin';
        s = pars.v1SpatialFilters;
        t = pars.v1TemporalFilters;
        dims = dims + [size(s,1), size(s,1), size(t,1)] - 1;

    case 'v1halfrect';
        s = pars.v1SpatialFilters;
        t = pars.v1TemporalFilters;
        dims = dims + [size(s,1), size(s,1), size(t,1)] - 1;

    case 'v1simple';
        s = pars.v1SpatialFilters;
        t = pars.v1TemporalFilters;
        dims = dims + [size(s,1), size(s,1), size(t,1)] - 1;

        if strcmp(pars.v1NormalizationType, 'untuned')
            dims = dims + [15 15 0];
        elseif strcmp(pars.v1NormalizationType, 'tuned1')
            s = pars.v1NormalizationSpatialFilter;
            t = pars.v1NormalizationTemporalFilter;
            dims = dims + [length(s)-1, length(s)-1, size(t,1)-1];
        elseif strcmp(pars.v1NormalizationType, 'tuned2')
            % not supported yet
        end

    case 'v1fullrect'
        s = pars.v1SpatialFilters;
        t = pars.v1TemporalFilters;
        dims = dims + [size(s,1), size(s,1), size(t,1)] - 1;
        
    case 'v1blur';
        s = length(pars.v1SpatialFilters) - 1;
        t = length(pars.v1TemporalFilters) - 1;
        dims = dims + [s s t];
        
        c = length(pars.v1ComplexFilter) - 1;
        dims = dims + [c c 0];

    case 'v1complex';
        s = length(pars.v1SpatialFilters) - 1;
        t = length(pars.v1TemporalFilters) - 1;
        dims = dims + [s s t];
        
        c = length(pars.v1ComplexFilter) - 1;
        dims = dims + [c c 0];
        
        if strcmp(pars.v1NormalizationType, 'global')
        else
            nx = length(pars.v1NormalizationSpatialFilter) - 1;
            nt = length(pars.v1NormalizationTemporalFilter) - 1;
            dims = dims + [nx nx nt];
        end


    case 'mtlin';
        s = length(pars.v1SpatialFilters) - 1;
        t = length(pars.v1TemporalFilters) - 1;
        dims = dims + [s s t];
        
        c = length(pars.v1ComplexFilter) - 1;
        dims = dims + [c c 0];
        
        if strcmp(pars.v1NormalizationType, 'global')
        else
            nx = length(pars.v1NormalizationSpatialFilter) - 1;
            nt = length(pars.v1NormalizationTemporalFilter) - 1;
            dims = dims + [nx nx nt];
        end

        
    case 'mtprepool';
        s = length(pars.v1SpatialFilters) - 1;
        t = length(pars.v1TemporalFilters) - 1;
        dims = dims + [s s t];
        
        c = length(pars.v1ComplexFilter) - 1;
        dims = dims + [c c 0];
        
        if strcmp(pars.v1NormalizationType, 'global')
        else
            nx = length(pars.v1NormalizationSpatialFilter) - 1;
            nt = length(pars.v1NormalizationTemporalFilter) - 1;
            dims = dims + [nx nx nt];
        end


        if pars.mtSpatialPoolingBeforeThreshold == 1
            s = length(pars.mtSpatialPoolingFilter) - 1;
            dims = dims + [s s 0];
        end
            
    case 'mthalfrect';
        s = length(pars.v1SpatialFilters) - 1;
        t = length(pars.v1TemporalFilters) - 1;
        dims = dims + [s s t];
        
        c = length(pars.v1ComplexFilter) - 1;
        dims = dims + [c c 0];
        
        if strcmp(pars.v1NormalizationType, 'global')
        else
            nx = length(pars.v1NormalizationSpatialFilter) - 1;
            nt = length(pars.v1NormalizationTemporalFilter) - 1;
            dims = dims + [nx nx nt];
        end


        if pars.mtSpatialPoolingBeforeThreshold == 1
            s = length(pars.mtSpatialPoolingFilter) - 1;
            dims = dims + [s s 0];
        end
            
    case 'mtpostpool';
        s = length(pars.v1SpatialFilters) - 1;
        t = length(pars.v1TemporalFilters) - 1;
        dims = dims + [s s t];
        
        c = length(pars.v1ComplexFilter) - 1;
        dims = dims + [c c 0];
        
        if strcmp(pars.v1NormalizationType, 'global')
        else
            nx = length(pars.v1NormalizationSpatialFilter) - 1;
            nt = length(pars.v1NormalizationTemporalFilter) - 1;
            dims = dims + [nx nx nt];
        end


        if pars.mtSpatialPoolingBeforeThreshold == 1
            s = length(pars.mtSpatialPoolingFilter) - 1;
            dims = dims + [s s 0];
        end
   
        if pars.mtSpatialPoolingBeforeThreshold == 0
            s = length(pars.mtSpatialPoolingFilter) - 1;
            dims = dims + [s s 0];
        end
        
    case 'mtpattern';
        s = length(pars.v1SpatialFilters) - 1;
        t = length(pars.v1TemporalFilters) - 1;
        dims = dims + [s s t];
        
        c = length(pars.v1ComplexFilter) - 1;
        dims = dims + [c c 0];
        
        if strcmp(pars.v1NormalizationType, 'global')
        else
            nx = length(pars.v1NormalizationSpatialFilter) - 1;
            nt = length(pars.v1NormalizationTemporalFilter) - 1;
            dims = dims + [nx nx nt];
        end

        if pars.mtSpatialPoolingBeforeThreshold == 1
            s = length(pars.mtSpatialPoolingFilter) - 1;
            dims = dims + [s s 0];
        end
   
        if pars.mtSpatialPoolingBeforeThreshold == 0
            s = length(pars.mtSpatialPoolingFilter) - 1;
            dims = dims + [s s 0];
        end
        
        if strcmp(pars.mtNormalizationType, 'global')
        else
            xSz = length(pars.mtNormalizationSpatialFilter) - 1;
            tSz = length(pars.mtNormalizationTemporalFilter) - 1;
            dims = dims + [xSz xSz tSz];
        end
        
    otherwise
        error([stageName, ' is not a recognized stage name.']);
end
