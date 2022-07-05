function [finalSettings,badIndex,quantized,perError,settings] = LinearToSettingsAcc(cal,linear)% [finalSettings,badIndex,quantized,perError,settings] = LinearToSettingsAcc(cal,linear)%% Convert from linear color space coordinates to device% setting coordinates.  This routine makes use of the% full basis function information to compensate for spectral% shifts in the device primaries with input settings.%% This depends on the standard calibration globals.%% 11/12/93   dhb      Wrote it.% 3/30/94	 dhb, jms Fixed logic bug in error computation.%					  Return finalSettings as best during iteration% 8/4/96     dhb      Update for stuff bag routines.% 8/21/97    dhb      Update for structures.% 3/10/98	 dhb	  Change nBasesOut to nPrimaryBases.% Algorithm parametersnIterations = 10;dampingFactor = 1.0;% Determine sizes[nLinear,nTargets] = size(linear);if (nTargets ~= 1)	error('Only handles one linear target at a time');endsettings = zeros(nLinear,nIterations);quantized = zeros(nLinear,nIterations);error = zeros(nLinear,nIterations);% Get basis informationnPrimaryBases = cal.nPrimaryBases;if (isempty(nPrimaryBases))	error('No nPrimaryBases field present in calibration structure');end% THINK ABOUT OUT OF GAMUT ISSUE.  THIS COMMENTED% OUT CODE WAS AN INITIAL STAB AT IT.%device = LinearToDevice(cal,linear);%[nDevice,null] = size(device);%gamut = DeviceToGamut(cal,device);%target = DeviceToLinear(cal,gamut);target = linear;aimfor = target;for i = 1:nIterations   device = LinearToDevice(cal,aimfor);  [gamut,badIndex] = DeviceToGamut(cal,device);  settings(:,i) = GamutToSettings(cal,gamut);  [tmpQuantized,deviceE] = SettingsToLinearAcc(cal,settings(:,i));	quantized(:,i) = tmpQuantized;  calError(:,i) = quantized(:,i) - aimfor;	perError(:,i) = quantized(:,i) - target;  aimfor = target - dampingFactor*calError(:,i);end% Find minimum error that was encountered and return those settingssummaryError = diag(perError'*perError);[null,minIndex] = min(summaryError);finalSettings = settings(:,minIndex);