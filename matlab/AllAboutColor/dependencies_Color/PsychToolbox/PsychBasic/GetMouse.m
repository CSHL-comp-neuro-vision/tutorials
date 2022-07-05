function GetMouse% [x,y,button] = GetMouse([windowPtrOrScreenNumber])% % Return the current (x,y) position of the cursor and the up/down state of% the mouse button. "button" is true (1) while pressed, and false (0) otherwise.%% The cursor position (x,y) is "local", i.e. relative to the origin of% the window or screen, if supplied. Otherwise it's "global", i.e. relative% to the origin of the main screen (the one with the menu bar).% % GetMouse also supports an obsolete usage:% xy = GetMouse([windowPtrOrScreenNumber])% where xy is a 1x2 vector containing the x, y coordinates.% % NOTE: If you use GetMouse to wait for clicks, then don't forget to wait% for the user to release the mouse button, ending the current click, before% you begin waiting for the next mouse press.% 4/27/96  dhb  Wrote this help file.% 5/12/96  dgp  Removed confusing comment about columns.%               Added query about coordinates.% 5/16/96  dhb  Modified MEX file to conform to above usage, answered%               query about coordinates.% 5/29/96  dhb  Flushing mouse button events added by dgp.% 8/23/96  dhb  Added support for windowInfo argument.% 2/24/97  dgp	Updated.% 2/24/97  dgp	Updated comments about flushing mouse button events.% 3/10/97  dgp	windowPtrOrScreenNumber% 3/23/97  dgp	deleted obsolete comment about flushing mouse button events.% 5/9/00   dgp  Added note about waiting for release before waiting for next click.