function [spd,S] = MeasMonSpd(window,settings,S,syncMode,bits,noMeterAvail,RADIUS)% [spd,S] = MeasMonSpd(window,settings,[S],[syncMode],[bits],[noMeterAvail],[RADIUS])%% Measure the Spd of a series of monitor settings.%% This routine is specific to go with CalibrateMon,% as it depends on the action of SetMon. %% If noMeterAvail is passed and set to 1, then the routine% returns random spectra.  This is useful for testing when% you don't have a meter.%% 10/26/93  dhb		Wrote it based on ccc code.% 11/12/93  dhb		Modified to use SetColor.% 8/11/94		dhb		Sync mode.% 8/15/94		dhb		Sync mode as argument, allow S to be [] for default.% 4/12/97   dhb   New toolbox compatibility, take window and bits args.% 8/26/97   dhb, pbe Add noMeterAvail arg.% 4/7/99    dhb   Add argument for radius board. Compact default arg code.% Check args and make sure window is passed right.usageStr = 'Usage: [spd,S] = MeasMonSpd(window,settings,[S],[syncMode],[bits],[noMeterAvail])';if (nargin < 2 | nargin > 7 | nargout > 2)	error(usageStr);endif (size(window,1) ~= 1 | size(window,2) ~= 1)	error(usageStr);end% Set defaultsdefaultS = [380 5 81];defaultSync = 0;defaultBits = 8;defaultNoMeterAvail = 0;% Check args and set defaultsif (nargin < 7 | isempty(RADIUS))	RADIUS = 0;endif (nargin < 6 | isempty(noMeterAvail))	noMeterAvail = defaultNoMeterAvail;endif (nargin < 5 | isempty(bits))	bits = defaultBits;endif (nargin < 4 | isempty(syncMode))	syncMode = defaultSync;endif (nargin < 3 | isempty(S))	S = defaultS;end% Radius default gammaif (RADIUS)	radiusGamma = SCREEN(window,'Gamma');end[null,nMeas] = size(settings);spd = zeros(S(3),nMeas);for i=1:nMeas	%fprintf(1,'Measurement %g ...',i);	if (RADIUS)		radiusGamma(2,:) = settings(:,i)';		SCREEN(window,'WaitBlanking');		SCREEN(window,'Gamma',radiusGamma,10);	else  	SetColor(window,1,settings(:,i),bits);	end	%fprintf(1,'set color ...');		% Set parameters for monitor measurements	if (~noMeterAvail)		if (syncMode ~= 0)			CMETER('SetParams',1,0,1);		else			CMETER('SetParams',0,0,1);		end	end	% Measure spectrum	if (~noMeterAvail)  	spd(:,i) = MeasSpd(S);	else		spd(:,i) = sum(settings(:,i))*ones(S(3),1);		WaitSecs(.1);	endend