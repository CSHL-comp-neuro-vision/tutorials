function directory=CalDataFolder% directory=CalDataFolder%% Get the path to the CalData folder.% Denis Pelli 7/25/96% Denis Pelli 2/28/98 change "CalDat" to "PsychCalData"name='PsychCalData';directory=FindFolder(name);if isempty(directory)	error(['Can''t find any ''' name ''' folder in the MATLAB path.']);endif size(directory,1)>1	for i=1:size(directory,1)		disp(['DUPLICATE: ''' deblank(directory(i,:)) '''']);	end	error(['Found more than one ''' name ''' folder in the MATLAB path.']);end