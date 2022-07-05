% MenuBarTest%% This little program demonstrates the flashing menu bar% problem when the new (Apple recommended) HideMenuBar is used, with% OS earlier than 9, with Rush.% % The symptom is that the menu bar region flashes during CLUT% animation.  The problem is also seen with ContrastThreshDemo, but it% is sufficiently complicated that it seemed worthwhile to have a more% direct test.% % This test program was used to determine the best default settings% for how to hide the menu bar (new vs old way). We believe that the% default settings (new in Mac OS 9 or later, old in older systems)% will be trouble free. However, if running this program reveals a% problem with the menu bar, you might want to try running it again% after using% SCREEN('Preference','UseNewHideMenuBar',0)% or% SCREEN('Preference','UseNewHideMenuBar',1)% to try the other setting.% 06/xx/00  dhb  Wrote it.% 07/7/00   dgp  Added dhb's email as a comment.% 07/09/00  dhb  Cleaned up, added clearer comments.% 07/10/00  dgp  Cosmetic.% 07/10/00  dhb  Add FORCE_USE_NEW% 07/10/00  dgp  Renamed FORCE_USE_NEW to useNewHideMenuBarclear allSHOW_NON_RUSH = 0;USE_EVAL_NOT_RUSH = 1;useNewHideMenuBar = []; % set this to 1, 0, or []. % [] requests default, which is 1 for Mac OS 9 or later, and 0 for earlier systems.% Setting it to 1 reveals a bug in Mac OS 8.6.% Set user's preference.oldUseNewHideMenuBar=SCREEN('Preference','UseNewHideMenuBar',useNewHideMenuBar);% Report current setting.comp=SCREEN('Computer');fprintf('%s\n',comp.system);fprintf('Preference UseNewHideMenuBar %d.\n',SCREEN('Preference','UseNewHideMenuBar')); % Set up some CLUTSnCluts = 80;theCluts = zeros(256,3,nCluts);theContrasts = sin(2*pi*(0:nCluts-1)/nCluts);for i = 1:nCluts	contrast = theContrasts(i);	if contrast==0		contrast=0.00001;	end	lowValDev = (1-contrast)/2;	highValDev = (1+contrast)/2;	clutEntries = round(linspace(lowValDev*255,highValDev*255,256)');	theCluts(:,:,i) = [clutEntries clutEntries clutEntries];end% Open screen and write in a simple image.[w,r] = SCREEN(0,'OpenWindow',128,[],8);SCREEN(w,'PutImage',(0:255)'*ones(1,128));% Clut movie without rush.  Simply contrasts reverses the% clut.  No problem.if SHOW_NON_RUSH	clutCounter = 1;	for j=1:nCluts		SCREEN(w,'SetClut',theCluts(:,:,clutCounter),0);		clutCounter=rem(clutCounter,nCluts)+1;		SCREEN(w,'WaitBlanking',2);	endend% Rush/eval the same movie.  Actually, eval() is sufficient to produce% the problem, a flashing menu bar region during CLUT animation.clutCounter = 1;s = ['for j=1:nCluts;' ...				'SCREEN(w,''SetClut'',theCluts(:,:,clutCounter),0);' ...				'clutCounter=rem(clutCounter,nCluts)+1;' ...				'SCREEN(w,''WaitBlanking'',2);' ...			'end;'];p=MaxPriority(w,'WaitBlanking','GetSecs');if USE_EVAL_NOT_RUSH	eval(s);else	Rush(s,p);endSCREEN(w,'Close');% RestoreSCREEN('Preference','UseNewHideMenuBar',oldUseNewHideMenuBar);