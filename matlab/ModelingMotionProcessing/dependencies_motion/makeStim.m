function stim = makeStim(xProfile,tProfile)

stim = zeros(length(tProfile),length(xProfile));

tProfile = round(tProfile);

for row=1:length(tProfile)
  shift = mod(tProfile(row)-1,length(xProfile))+1;
  stim(row,:) = circularShift(xProfile,shift,0);
end
