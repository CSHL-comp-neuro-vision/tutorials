function bounds=TextBounds(w,text)% bounds=TextBounds(window,string)%% Returns the smallest enclosing rect for the drawn text, relative to% the current location. This bound is based on the actual pixels% drawn, so it incorporates effects of text smoothing, etc. "text"% may be a cell array or matrix of 1 or more strings. The strings are% drawn one on top of another, at the same initial position, before% the bounds are calculated. This returns the smallest box that will% contain all the strings. The prior contents of the scratch window% are lost. Usually it should be an offscreen window, so the user% won't see it. The scratch window should be at least twice as wide% as the text, because the text is drawn starting at the left-right% center of the window to cope with uncertainties about text% direction (e.g. Hebrew) and some unusual characters that extend% greatly to the left of their nominal starting point. If you only% know your nominal text size and number of characters, you might do% this to create your scratch window:% w=SCREEN(0,'offscreenwindow',[],[0 0 3*size*chars 3*size]);% % Also see TextCenteredBounds.% 9/1/98 dgp wrote it.% 3/19/00 dgp debugged it.SCREEN(w,'FillRect',0);r=SCREEN(w,'Rect');x0=(r(RectLeft)+r(RectRight))/2;y0=round((r(RectTop)+2*r(RectBottom))/3);for i=1:size(text,1)	string=char(text(i,:));	SCREEN(w,'DrawText',string,x0,y0,255);endSCREEN(w,'DrawText','',x0,y0);image1=SCREEN(w,'GetImage');[y,x]=find(image1);if isempty(y) | isempty(x)	bounds=[0 0 0 0];else		bounds=SetRect(min(x)-1,min(y)-1,max(x),max(y));	bounds=OffsetRect(bounds,-x0,-y0);end