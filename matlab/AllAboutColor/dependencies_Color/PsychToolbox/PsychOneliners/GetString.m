function string = GetString% string = GetString% % Get a string typed at the keyboard. Entry is terminated by % <return> or <enter>.%% Useful for i/o in a SCREEN window. Typed keys are not echoed.%% See also: GetEchoString, GetNumber, GetEchoNumber% 12/7/95 dhb	Wrote GetNumber in response to query from Tina Beard.% 12/8/95 dhb	Add delete functionality.% 2/4/97  dhb	Fixed bug.  Can now hit delete more than once.% 2/5/97  dhb	Accept <enter> as well as <cr>.%         dhb	Allow string return as well.  % 3/15/97 dgp Created GetString based on dhb's GetNumber.% 3/31/97 dhb Fix bug arising from new initialization.% 2/28/98 dgp Use GetChar instead of obsolete GetKey. Use SWITCH and LENGTH.% 3/27/98 dhb Fix bug from 2/28/98, put abs around char in switch.string='';while 1	% Loop until <return> or <enter>	char=GetChar;	switch(abs(char))		case {13,3},	% <return> or <enter>			break;		case 8,			% <delete>			if length(string)>0				string=string(1:length(string)-1);			end		otherwise,			string=[string char];	endend