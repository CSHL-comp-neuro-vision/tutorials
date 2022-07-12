%% Parfor example
% by Dan Birman (2022)
%
% There is a lot of overhead to using a parfor loop in matlab, so it's best
% to have each step of the parfor loop take a few seconds or more to
% calculate. If your steps are too fast, you can "chunk" your data to make
% each parfor step go more slowly.

% define a matrix of data with 495 variables (e.g. voxels), each of which
% has some amount of data attached (timepoints?)
data = ones(495,10000);

%% Default parfor
parfor vi = 1:size(data,1)
    loopData = data(vi,:);
    % do stuff on loopData
    data(vi,:) = loopData;
end

%% Chunked parfor
% lots of data setup
chunkSize = 100;
% figure out which chunk will use which data
chunkStarts = 1:(chunkSize-1):size(data,1);
% generate the chunk indexes across the first dimension of data
chunkIdxs = [];
for ci = 1:(length(chunkStarts)-1)
    chunkIdxs(ci,:) = [chunkStarts(ci) chunkStarts(ci+1)];
end

% chunk the data into a new array that is nchunks size in the first dim
chunkData = zeros(size(chunkIdxs,1),chunkSize,size(data,2));
for ci = 1:size(chunkIdxs,1)
    chunk = data(chunkIdxs(ci,1):chunkIdxs(ci,2),:);
    
    % we have to do some funk here because the final chunk might not be
    % full size
    chunkData(ci,:,:) = [chunk ; zeros(max(0,chunkSize-size(chunk,1)-1),size(chunkData,3))];
end

% now run the parfor loop acros schunks
parfor ci = 1:size(chunkData,1)
    % pull the chunk data
    localData = chunkData(ci,:,:);
    
    % loop over the chunk and do your thing, same as above
    for vi = 1:size(localData,1)
        loopData = localData(vi,:);
        % do stuff on loopData
        localData(vi,:) = loopData;
    end
    
    % save the output
    chunkData(ci,:,:) = localData;
end