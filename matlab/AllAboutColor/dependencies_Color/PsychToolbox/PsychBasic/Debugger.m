function Debugger(string)% Debugger([string])% % Pass control to whatever debugger is installed in your Macintosh. If no% debugger is installed, then Debugger.mex issues an error message. (Until % 12/99, calling Debugger.mex when no debugger was installed would crash % with System error 12.) See below to download Apple's free MacsBug. Once % you're in the debugger, most debuggers will accept the go command "g", % followed by Return, to resume your application.% % The optional "string" (up to 255 characters) is passed to the debugger,% to remind you how you got there, or to issue commands to the debugger.% Precede each MacsBug command by a semicolon. E.g.% % Debugger('Hello Hades!')% displays a hello message. Type G Return when you're ready to continue.% % Debugger('Hello Hades!;g')% briefly displays the hello and returns immediately.% % Debugger(';hc all')% checks all the heaps. Type G Return when you're ready to continue.% % for i=1:10% 	Debugger(['i=' num2str(i) ';hc;g;']);% end% will do ten heap checks, each labeled for your review when you next% examine the MacsBug display.% % For further information on MacsBug, invoke MacsBug (e.g. by typing% Debugger in MATLAB) and then type question mark "?" followed by Return.% You can get help on any one topic by typing ? followed by the topic% name. When you're satisfied, type "g", followed by Return, to resume% MATLAB. (MacsBug ignores upper/lower case in most instances.)% % MacsBug has �cw� and �cwp� macros to take you to the Metrowerks% CodeWarrior debugger in 68K and PowerPC code, respectively.% % When your computer crashes, and you end up in MacsBug, if you want to% write to the author of the program that you think caused the crash, it's% very helpful to document the machine's state and recent history (i.e.% routines on the stack). The MacsBug command "StdLog" saves this% diagnostic information as a file called "StdLog" on your desktop, which% you can mail to the programmer.% % MacsBug is free. Select the following line (triple click) and hit enter:% web http://developer.apple.com/tools/debuggers/MacsBug/% % More MacsBug-related links:% web http://vision.nyu.edu/Tips/RecSoftware.html#MacsBug% % See LogTrapsTest.% % xx/xx/xx  dhb  Wrote it.% 5/28/96   dgp  Added string passing and expanded help.% 1/25/97   dgp  Updated MacsBug link address. Mention StdLog.% 1/26/97   dgp  Added another example; mention cw,cwp.% 7/20/97   dgp  Added link to MacsBug info in VideoToolbox site.% 8/2/97	dgp  Updated url.% 3/16/98   dgp  Update.% 12/23/99	dgp  Updated link to apple's web site.