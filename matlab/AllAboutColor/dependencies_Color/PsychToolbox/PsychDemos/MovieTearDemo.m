function MovieTearDemo% MovieTearDemo%% Drifting vertical edge demonstrates "tearing".% % TEARING: You may have seen "tearing" in the movie. The tear is a% horizontal break, perhaps a third of the way down from the top of the% screen, jiggling up and down as the movie played. This tear is normal% and does not indicate anything wrong with your computer or the display% software. Tearing can occur in any display that isn't double buffered,% and hardly any of the drivers for Macintosh video cards support double% buffering. Even so, it's easy to avoid tearing once you understand why% it happens. % % The tear is a discontinuity in time. The movie above the tear is delayed% by one frame relative to what's displayed below the tear. This is% invisible if the movie is static, but produces a glaring break in a% moving object. Your video card has enough video memory to store one% frame. The reading of video memory is an autonomous free-running process% that drives the intensity and position of the monitor's video beam in a% raster scan that paints the screen from top to bottom. The tear is the% result of the computer writing to the same video memory that the video% beam is reading. The computer is racing with the video beam from the top% to the bottom of your screen. The beam has a head start, but your% computer is slightly faster. The tear occurs when your computer passes% the beam. Below the tear you see the contents of the new frame, just% written by your computer into the frame store. Above the tear you are% seeing the contents of the previous frame. % % You can eliminate the tear by pushing it up off the top of your movie.% Anything that makes your computer proceed more quickly down the screen% (e.g. making the movie narrower) will move the tear higher, as will% moving your movie's window lower, i.e. giving your computer a head start% in the race from top to bottom of the screen.% % JIGGLING: The relatively stable vertical position of the tear is a% consequence of the synchrony of the computer and monitor; if they% weren't synchronized the tear would occur at a different height on each% frame (which would make it much harder to get rid of). The jiggling up% and down of the tear reveals the computer's varying latency, due to Mac% OS and device driver interrupts. This allows you to directly observe the% effect of using RUSH to minimize interruptions.% % DOUBLE BUFFERING: The hardware method of eliminating tearing is double% buffering; you use two pages of video memory and alternately display% from one while writing the other. The Mac OS defines driver-level calls% to support double buffering (multiple pages) and most video cards have% suitable chips and enough video memory to provide multiple pages, but% hardly any video-card manufacturers have bothered to program their video% drivers to support more than one page of video memory. ScreenTest% reports how many video pages you have. We could add MATLAB-level calls% to control double buffering. Please send me a note if you have a video% device with more than one video page and you'd like to try double% buffering. denis@psych.nyu.edu% 2/4/98	dgp	 Wrote it, based on MovieDemo. Show sliding vertical line %							 instead of growing disk.% 2/16/98 dgp	 Use edge instead of line. Rewrite explanation.% 3/10/98	dhb	 Load DefaultClut.% 7/14/98	dgp	 ClutDefault instead of DefaultClut.% 7/17/98	dgp   Using enhanced Rush, use easy-to-read cell array for string.% Open an onscreen windowframes=128;screenNumber=0;pixelSize=8;rect=SCREEN(screenNumber,'Rect');rect(RectRight)=128;window=SCREEN(screenNumber,'OpenWindow',0,rect,pixelSize);SCREEN(window,'SetClut',ClutDefault(window));% Move MATLAB Command window to the right of oursm=SCREEN('GetMATLABWindow');r=AdjoinRect(SCREEN(m,'GlobalRect'),SCREEN(window,'GlobalRect'),RectRight);r=OffsetRect(r,7,0);SCREEN(m,'GlobalRect',r);% Make a movie by drawing an edge into 128 offscreen windows.fprintf('Computing the movie... ');r=OffsetRect(rect,1-RectWidth(rect),0);step=max(1,round((RectWidth(rect)-RectWidth(r))/(frames-1)));imagePtr=zeros(frames,1);for i=1:frames	imagePtr(i)=SCREEN(window,'OpenOffscreenWindow',0);	SCREEN(imagePtr(i),'FillRect',255,OffsetRect(r,(i-1)*step,0));endfprintf('Done.\n');HideCursor;% Show the movie by looping on CopyWindow.% The movie is shown forwards, then backwards.for i=[1:frames frames:-1:1]	SCREEN(window,'WaitBlanking');	SCREEN('CopyWindow',imagePtr(i),window);end% Show the movie again, now using Rush to minimize interruptions.loop={	'for i=[1:frames frames:-1:1];'		'SCREEN(window,''WaitBlanking'');'		'SCREEN(''CopyWindow'',imagePtr(i),window);'	'end;'};priorityLevel=MaxPriority(screenNumber,'WaitBlanking');SCREEN('Screens'); % Make sure all Rushed functions are in memory.Rush(loop,priorityLevel);% Close all the on- and off-screen windowsSCREEN('CloseAll');ShowCursor;fprintf('The movie was shown twice. The second showing used RUSH to minimize\n');fprintf('interruptions. See RUSH.\n\n')% frameRate[frame0,s0]=screen(screenNumber,'PeekBlanking');screen(screenNumber,'WaitBlanking',2);[frame,s]=screen(screenNumber,'PeekBlanking');frameRate=(frame-frame0)/(s-s0);% secsPerLiner=SCREEN(screenNumber,'Rect');secsPerLine=((1/frameRate)-0.001)/RectHeight(r);	% allow for 1 ms of blankingfprintf('TEARING: You may have seen "tearing" in the movie. The tear is a\n');fprintf('horizontal break, perhaps a third of the way down from the top of the\n');fprintf('screen, jiggling up and down as the movie played. This tear is normal\n');fprintf('and does not indicate anything wrong with your computer or the display\n');fprintf('software. Tearing can occur in any display that isn''t double buffered,\n');fprintf('and hardly any of the drivers for Macintosh video cards support double\n');fprintf('buffering. Even so, it''s easy to avoid tearing once you understand why\n');fprintf('it happens. For an explanation:\n');fprintf('help MovieTearDemo	% select this line and hit the <enter> key.\n');fprintf('\n');fprintf('JIGGLING: The vertical jiggling of the tear tracks your computer''s\n');fprintf('synchronization to the display. Each pixel of jiggle up and down\n');fprintf('corresponds to %.0f microseconds change in the synchrony of computer and\n',secsPerLine*1e6);fprintf('display. RUSH improves synchrony, reducing the jiggle. See RUSH.\n');