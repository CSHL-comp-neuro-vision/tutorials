% BytesTest.m%% Demonstrate the Bytes mex file. %% dgp 1/20/96  Wrote it% dgp 3/16/97  Updatedfree=system_dependent(1003)			% undocumented built-in MATLAB functionfree=Bytes							% Apple's FreeMem(). Innocuous.free=Bytes('Free')					% "          "          "tempFree=Bytes('TempFree')			% Apple's TempFreeMem(). Innocuous.maxBlock=Bytes('MaxBlock')			% Apple's MaxBlock(). Innocuous.stackSpace=Bytes('StackSpace')		% Apple's StackSpace(). Innocuous.disp('WARNING: The following two calls may take a while; they move memory.')maxMem=Bytes('MaxMem')				% Apple's MaxMem(). CAUTION: moves memory!tempMaxMem=Bytes('TempMaxMem')		% Apple's TempMaxMem(). CAUTION: moves memory!