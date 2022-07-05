function wls = MakeItWls(S)% wls = MakeItWls(S)%% If argument is a [start delta n] description, it is% expanded to an actual list of wavelengths.%% If passed length is not a [start delta n] description,% it is checked to see if it is a column vector with evenly% spaced entries.% Check if we are a [start delta n][m,n] = size(S);if (m == 1 & n == 3)  wls = SToWls(S);% If not, assume its already a list of wavelengths.  But check.else  % As a side effect, the call to WlsToS checks that the passed  % description is a valid list of wavelengths  wls = S;  WlsToS(wls);end