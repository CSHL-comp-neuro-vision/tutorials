function wls = SToWls(S)% wls = SToWls(S)%% Expand a [start delta n] description to an actual% list of wavelengths.% Check validity[m,n] = size(S);if (m ~= 1 | n ~= 3)  error('Passed list is not a [start delta n] description');endif (S(1) <= 0 | S(2) == 0 | S(3) <=0)  error('Passed list is not a [start delta n] description');end% Expand awaywls = (S(1):S(2):S(1)+(S(3)-1)*S(2))';