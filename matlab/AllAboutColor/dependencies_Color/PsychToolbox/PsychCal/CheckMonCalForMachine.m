function isOK = CheckMonCalForMachine(cal,whichScreen)% isOK = CheckMonCalForMachine(cal,whichScreen)%% Check that the information in the passed calibration% file is consisent with the current computer, driver, etc.%% 8/26/97  dhb, pbe  Wrote it.% 9/18/97  dhb, pbe  Change name sign convention on return.% 3/10/98	 dhb			 Change name to CheckMonCalForMachine.% 7/3/98   dhb, pbe  Change for cal.describe format.% Fill in descriptive information[computerName,owner,system,processor,cache,fpu,Hz,busHz,vm,pci]=SCREEN('Computer');checkComputer = sprintf('%s''s %s, %s, %s',owner,computerName,system,vm);checkScreen = whichScreen;[card,driver,driverVersion,slot]=SCREEN(whichScreen,'VideoCard');checkDriver = sprintf('%s (%s version %g) in slot %s',card,driver,driverVersion,slot);checkDacsize = SCREEN(whichScreen,'Preference','ClutDacSize');isOK = 1;if (~strcmp(cal.describe.computer,checkComputer))	%fprintf(1,'CheckCalForMachine:\n');	%fprintf(1,'\tCurrent computer:     %s',checkComputer);	%fprintf(1,'\tCalibration computer: %s',cal.describe.computer);	isOK = 0;endif (cal.describe.whichScreen ~= checkScreen)	%fprintf(1,'CheckCalForMachine:\n');	%fprintf(1,'\tCurrent screen:       %g\n',checkScreen);	%fprintf(1,'\tCalibration computer: %g\n',cal.describe.whichScreen);	isOK = 0;end	if (~strcmp(cal.describe.driver,checkDriver))	%fprintf(1,'CheckCalForMachine:\n');	%fprintf(1,'\tCurrent driver:     %s\n',checkDriver);	%fprintf(1,'\tCalibration driver: %s\n',cal.driver);	isOK = 0;endif (cal.describe.dacsize ~= checkDacsize)	%fprintf(1,'CheckCalForMachine:\n');	%fprintf(1,'\tCurrent DAC size:       %g\n',checkDacsize);	%fprintf(1,'\tCalibration DAC size:   %g\n',cal.describe.dacsize);	isOK = 0;end