function LogTraps% LogTraps% Records all the Apple traps executed while running a bit of MATLAB% code. This includes traps in any MEX files that are executed.% The log is saved in a text file "TrapLog" on your desktop.% You must have MacsBug installed. LogTraps is useful for 68K code.% In PowerPC code, all traps are compiled as subroutine calls, which% are missed by MacsBug's trap-recording mechanism.%% You can create a 68K-only version of MATLAB, from a fat original,% by using ResEdit to remove the "cnfg" resource. Do this to a copy;% don't modify your original application. To run it on a PowerMac% you'll need a floating point emulator, e.g. SoftwareFPU.%% At present the bit of MATLAB code is 3 iterations of a FOR loop,% for i=1:n% 	GetSecs(t,i);% end% but you could edit the file to log traps in any MATLAB code.% 4/12/97 dgp wrote it% 4/13/97 dgp expanded comments% 4/15/97 dgp expanded commentsfprintf('Don''t be alarmed when I intentionally drop into MacsBug\n');fprintf('and noisily log a file to disk.\n');fprintf('You MUST have MacsBug installed to run this program!!\n');if ~input('Okay to continue (1=yes,0=no)?');	break;endfprintf('\n');% Describe the computer[computerName,owner,system]=SCREEN('Computer');id=['MATLAB ' version ', ' system ', ' owner '''s ' computerName ];fprintf('%s\n',id);% Turn off things that might slow us down.fsWasOn=FileShare(-3,1);adWasOn=AfterDark(0);if FileShare ~= -3	fprintf('FileSharing is on, which may slow things down.\n');endfprintf('\n');% Load everything before we starti=0;n=2;t=zeros(1,n);GetSecs;% Log to diskdebugger([id ';atc;log TrapLog;set suspendprompt on;stat;atta ,0;g']); for i=1:n	%debugger('next i;set suspendprompt on;g')	GetSecs(t,i);	%bitand(i,i);enddebugger(';log;atc;g')% RestoreFileShare(fsWasOn,0);fprintf('Done. Results are in the "TrapLog" file on your desktop.\n');