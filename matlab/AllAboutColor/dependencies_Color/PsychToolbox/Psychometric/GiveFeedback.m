function GiveFeedback(correct)% GiveFeedback(correct)%% Give auditory feedback about correctness of respose.%% 3/5/97  	dhb  Wrote it% 1/25/00 	emw  Added platform conditionals% 3/8/2000	dgp  Fixed platform conditionals% 4/14/00   dhb  Fix call to SND for windows.if correct   if strcmp(computer,'PCWIN')      SND('Play',sin(1:1000));   else      SndPlay(sin((1:200)));   end   else   if strcmp(computer,'PCWIN')      SND('Play',sin((1:1000)));      WaitSecs(.2);      SND('Play',sin((1:1000)));   else      SndPlay(sin((1:200)));      WaitSecs(0.1);      SndPlay(sin((1:200)));   endend