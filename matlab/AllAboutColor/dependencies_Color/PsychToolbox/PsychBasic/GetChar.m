function character = GetChar% character = GetChar% % Wait (if necessary) for a typed character and return it.% % Command-Period always causes an immediate exit.% % GetChar and CharAvail are character-oriented (and slow), whereas% KbCheck and KbWait are keypress-oriented (and fast). If only a meta% key (like <option> or <shift>) was hit, KbCheck will return true,% because a key was pressed, but CharAvail will return false, because% no character was generated. See KbCheck.% % CharAvail and GetChar use the Event Manager to retrieve the% character generated, not the raw key press(es) per se. If the user% presses "a", GetChar returns 'a', but if the user presses option-e% followed by "a", this selects an accented a, "�", which is treated% by GetChar as a single character, even though it took the user three% keypresses (counting the option key) to produce it.% % CharAvail and GetChar call the Event Manager, which allows the% system to get control. Sometimes CharAvail will take tens of% milliseconds to return, so don't use CharAvail in real-time loops.% And there can be some delay between when the key is pressed and when% CharAvail or GetChar detects it. If precise timing of the keypress% is important, use KbCheck or KbWait.% % WARNING: When BACKGROUNDING is enabled, MATLAB removes all% characters from the event queue before executing each MATLAB% statement, so CharAvail and EventAvail('keyDown') always report 0.% So turn off BACKGROUNDING:% % SCREEN('Preference','Backgrounding',0); % Until MATLAB 5.2.1, this call required a disk access, which is slow.% % On Windows, the BACKGROUNDING stuff is not applicable.%% See also: KbCheck, KbWait, CharAvail, KbDemo, EventAail, SCREEN Preference Backgrounding.% 5/7/96  dgp	Wrote this help file.% 1/22/97 dhb	Added comment and pointer to TIMER routines.% 3/6/97  dhb	References to KbWait, KbCheck.% 7/23/97 dgp	It's a character not a keypress.% 8/2/97  dgp	Explain difference between key and character. See KbCheck.% 2/7/98  dgp	Streamlined. Eliminated call to GetKey, since it's now GetChar.mex.% 3/24/98 dgp	Explain backgrounding and meta keys. Don't mention obsolete GetKey and KbHit.% 3/15/99 xmz	Put in some comment for Windows version.% 3/19/99 dgp	Update explanation of backgrounding. % 3/28/99 dgp	Show how to turn off backgrounding. 