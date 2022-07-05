function [ret,x0,str,ts,xts]=simulink_modules(t,x,u,flag);%SIMULINK_MODULES	is the M-file description of the SIMULINK system named SIMULINK_MODULES.%	The block-diagram can be displayed by typing: SIMULINK_MODULES.%%	SYS=SIMULINK_MODULES(T,X,U,FLAG) returns depending on FLAG certain%	system values given time point, T, current state vector, X,%	and input vector, U.%	FLAG is used to indicate the type of output to be returned in SYS.%%	Setting FLAG=1 causes SIMULINK_MODULES to return state derivatives, FLAG=2%	discrete states, FLAG=3 system outputs and FLAG=4 next sample%	time. For more information and other options see SFUNC.%%	Calling SIMULINK_MODULES with a FLAG of zero:%	[SIZES]=SIMULINK_MODULES([],[],[],0),  returns a vector, SIZES, which%	contains the sizes of the state vector and other parameters.%		SIZES(1) number of states%		SIZES(2) number of discrete states%		SIZES(3) number of outputs%		SIZES(4) number of inputs%		SIZES(5) number of roots (currently unsupported)%		SIZES(6) direct feedthrough flag%		SIZES(7) number of sample times%%	For the definition of other parameters in SIZES, see SFUNC.%	See also, TRIM, LINMOD, LINSIM, EULER, RK23, RK45, ADAMS, GEAR.% Note: This M-file is only used for saving graphical information;%       after the model is loaded into memory an internal model%       representation is used.% the system will take on the name of this mfile:sys = mfilename;new_system(sys)simver(1.3)if (0 == (nargin + nargout))     set_param(sys,'Location',[678,71,1009,473])     open_system(sys)end;set_param(sys,'algorithm',     'RK-45')set_param(sys,'Start time',    '0.0')set_param(sys,'Stop time',     '999999')set_param(sys,'Min step size', '0.0001')set_param(sys,'Max step size', '10')set_param(sys,'Relative error','1e-3')set_param(sys,'Return vars',   '')add_block('built-in/From Workspace',[sys,'/','Trigger'])set_param([sys,'/','Trigger'],...		'matl_expr','trigspec',...		'position',[35,329,75,351])add_block('built-in/From Workspace',[sys,'/',['Vertical',13,'Specification']])set_param([sys,'/',['Vertical',13,'Specification']],...		'orientation',3,...		'matl_expr','vspec',...		'position',[105,330,145,355])add_block('built-in/From Workspace',[sys,'/',['Horizontal',13,'Specification ']])set_param([sys,'/',['Horizontal',13,'Specification ']],...		'orientation',1,...		'matl_expr','hspec',...		'Mask Display','',...		'position',[220,325,260,355])%     Subsystem  ['Absolute-',13,'Time',13,'Specifier'].new_system([sys,'/',['Absolute-',13,'Time',13,'Specifier']])set_param([sys,'/',['Absolute-',13,'Time',13,'Specifier']],'Location',[433,163,667,398])add_block('built-in/From Workspace',[sys,'/',['Absolute-',13,'Time',13,'Specifier/Vertical',13,'Specification']])set_param([sys,'/',['Absolute-',13,'Time',13,'Specifier/Vertical',13,'Specification']],...		'matl_expr','vspec',...		'position',[50,167,90,193])add_block('built-in/From Workspace',[sys,'/',['Absolute-',13,'Time',13,'Specifier/Trigger']])set_param([sys,'/',['Absolute-',13,'Time',13,'Specifier/Trigger']],...		'matl_expr','trigspec',...		'position',[50,99,90,121])add_block('built-in/From Workspace',[sys,'/',['Absolute-',13,'Time',13,'Specifier/Horizontal',13,'Specification ']])set_param([sys,'/',['Absolute-',13,'Time',13,'Specifier/Horizontal',13,'Specification ']],...		'matl_expr','hspec',...		'Mask Display','',...		'position',[50,25,90,55])add_block('built-in/Mux',[sys,'/',['Absolute-',13,'Time',13,'Specifier/Mux']])set_param([sys,'/',['Absolute-',13,'Time',13,'Specifier/Mux']],...		'inputs','3',...		'position',[135,94,165,126])add_block('built-in/Outport',[sys,'/',['Absolute-',13,'Time',13,'Specifier/out_1']])set_param([sys,'/',['Absolute-',13,'Time',13,'Specifier/out_1']],...		'position',[185,100,205,120])add_line([sys,'/',['Absolute-',13,'Time',13,'Specifier']],[170,110;180,110])add_line([sys,'/',['Absolute-',13,'Time',13,'Specifier']],[95,110;130,110])add_line([sys,'/',['Absolute-',13,'Time',13,'Specifier']],[95,40;105,40;105,100;130,100])add_line([sys,'/',['Absolute-',13,'Time',13,'Specifier']],[95,180;105,180;105,120;130,120])%     Finished composite block ['Absolute-',13,'Time',13,'Specifier'].set_param([sys,'/',['Absolute-',13,'Time',13,'Specifier']],...		'position',[115,237,140,273])%     Subsystem  ['Relative-',13,'Time',13,'Specifier'].new_system([sys,'/',['Relative-',13,'Time',13,'Specifier']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier']],'Location',[321,301,788,790])add_block('built-in/MATLAB Fcn',[sys,'/',['Relative-',13,'Time',13,'Specifier/movementinc']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/movementinc']],...		'orientation',1,...		'MATLAB Fcn','movementinc(u(1))',...		'Output Width','1',...		'position',[47,265,93,295])%     Subsystem  ['Relative-',13,'Time',13,'Specifier/Timer'].new_system([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']],'Location',[638,313,1006,490])add_block('built-in/Note',[sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Outputs time from previous reset']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Outputs time from previous reset']],...		'position',[190,10,195,15])add_block('built-in/Sum',[sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Sum6']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Sum6']],...		'inputs','+-',...		'position',[250,74,265,101])add_block('built-in/Outport',[sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/time']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/time']],...		'position',[330,80,350,100])add_block('built-in/Inport',[sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/reset']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/reset']],...		'position',[20,65,40,85])add_block('built-in/Clock',[sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Clock1']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Clock1']],...		'Mask Display','',...		'position',[45,30,65,50])add_block('built-in/Memory',[sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Memory']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Memory']],...		'position',[95,90,130,120])add_block('built-in/Switch',[sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Switch']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer/Switch']],...		'Threshold','0.001',...		'position',[180,79,210,111])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']],[215,95;245,95])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']],[270,90;325,90])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']],[70,40;150,40;220,80;245,80])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']],[70,40;120,40;175,85])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']],[45,75;130,75;175,95])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']],[135,105;175,105])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']],[215,95;215,146;80,146;90,105])%     Finished composite block ['Relative-',13,'Time',13,'Specifier/Timer'].set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Timer']],...		'position',[45,192,65,218])add_block('built-in/Inport',[sys,'/',['Relative-',13,'Time',13,'Specifier/in_1']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/in_1']],...		'position',[10,195,30,215])add_block('built-in/Note',[sys,'/',['Relative-',13,'Time',13,'Specifier/Generates continuous movement specification and trigger signals in relative time.',13,'NOTE: must be used in conjunction with Movement Detector module']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Generates continuous movement specification and trigger signals in relative time.',13,'NOTE: must be used in conjunction with Movement Detector module']],...		'position',[260,20,265,25])add_block('built-in/MATLAB Fcn',[sys,'/',['Relative-',13,'Time',13,'Specifier/bighspec reader']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/bighspec reader']],...		'MATLAB Fcn','spec1(u(1))',...		'Output Width','1',...		'position',[200,61,245,99])add_block('built-in/Switch',[sys,'/',['Relative-',13,'Time',13,'Specifier/Switch']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Switch']],...		'position',[280,109,310,141])add_block('built-in/Memory',[sys,'/',['Relative-',13,'Time',13,'Specifier/Memory']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Memory']],...		'position',[195,120,230,150])add_block('built-in/Memory',[sys,'/',['Relative-',13,'Time',13,'Specifier/Memory1']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Memory1']],...		'position',[195,245,230,275])add_block('built-in/Switch',[sys,'/',['Relative-',13,'Time',13,'Specifier/Switch1']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Switch1']],...		'position',[280,234,310,266])add_block('built-in/MATLAB Fcn',[sys,'/',['Relative-',13,'Time',13,'Specifier/bigtspec reader']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/bigtspec reader']],...		'MATLAB Fcn','spec2(u(1))',...		'Output Width','1',...		'position',[200,186,245,224])add_block('built-in/Memory',[sys,'/',['Relative-',13,'Time',13,'Specifier/Memory2']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Memory2']],...		'position',[195,375,230,405])add_block('built-in/Switch',[sys,'/',['Relative-',13,'Time',13,'Specifier/Switch2']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Switch2']],...		'position',[280,364,310,396])add_block('built-in/MATLAB Fcn',[sys,'/',['Relative-',13,'Time',13,'Specifier/bigvspec reader']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/bigvspec reader']],...		'MATLAB Fcn','spec3(u(1))',...		'Output Width','1',...		'position',[200,316,245,354])add_block('built-in/Mux',[sys,'/',['Relative-',13,'Time',13,'Specifier/Mux']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/Mux']],...		'inputs','3',...		'position',[370,234,400,266])add_block('built-in/Outport',[sys,'/',['Relative-',13,'Time',13,'Specifier/out_1']])set_param([sys,'/',['Relative-',13,'Time',13,'Specifier/out_1']],...		'position',[420,240,440,260])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[70,205;70,260])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[70,205;160,205;160,335;195,335])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[70,205;160,205;160,80;195,80])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[70,205;195,205])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[35,205;40,205])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[315,380;325,380;325,430;175,430;175,390;190,390])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[250,335;260,335;260,380;275,380])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[250,335;260,335;260,370;275,370])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[235,390;275,390])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[235,260;275,260])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[315,250;315,300;180,300;190,260])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[250,205;260,205;260,240;275,240])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[250,205;260,205;260,250;275,250])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[250,80;260,80;260,125;275,125])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[250,80;260,80;260,115;275,115])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[315,125;315,175;180,175;190,135])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[235,135;275,135])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[405,250;415,250])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[315,125;335,125;335,240;365,240])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[315,250;365,250])add_line([sys,'/',['Relative-',13,'Time',13,'Specifier']],[315,380;335,380;335,260;365,260])%     Finished composite block ['Relative-',13,'Time',13,'Specifier'].set_param([sys,'/',['Relative-',13,'Time',13,'Specifier']],...		'position',[40,237,70,273])%     Subsystem  ['Pulse-Reset',13,'Timer'].new_system([sys,'/',['Pulse-Reset',13,'Timer']])set_param([sys,'/',['Pulse-Reset',13,'Timer']],'Location',[98,593,466,770])add_block('built-in/Switch',[sys,'/',['Pulse-Reset',13,'Timer/Switch']])set_param([sys,'/',['Pulse-Reset',13,'Timer/Switch']],...		'Threshold','0.001',...		'position',[180,79,210,111])add_block('built-in/Memory',[sys,'/',['Pulse-Reset',13,'Timer/Memory']])set_param([sys,'/',['Pulse-Reset',13,'Timer/Memory']],...		'position',[95,90,130,120])add_block('built-in/Clock',[sys,'/',['Pulse-Reset',13,'Timer/Clock1']])set_param([sys,'/',['Pulse-Reset',13,'Timer/Clock1']],...		'position',[45,30,65,50])add_block('built-in/Inport',[sys,'/',['Pulse-Reset',13,'Timer/reset']])set_param([sys,'/',['Pulse-Reset',13,'Timer/reset']],...		'position',[20,65,40,85])add_block('built-in/Outport',[sys,'/',['Pulse-Reset',13,'Timer/time']])set_param([sys,'/',['Pulse-Reset',13,'Timer/time']],...		'position',[330,80,350,100])add_block('built-in/Sum',[sys,'/',['Pulse-Reset',13,'Timer/Sum6']])set_param([sys,'/',['Pulse-Reset',13,'Timer/Sum6']],...		'inputs','+-',...		'position',[250,74,265,101])add_block('built-in/Note',[sys,'/',['Pulse-Reset',13,'Timer/Outputs elapsed time since previous reset']])set_param([sys,'/',['Pulse-Reset',13,'Timer/Outputs elapsed time since previous reset']],...		'position',[190,10,195,15])add_line([sys,'/',['Pulse-Reset',13,'Timer']],[215,95;215,146;80,146;90,105])add_line([sys,'/',['Pulse-Reset',13,'Timer']],[135,105;175,105])add_line([sys,'/',['Pulse-Reset',13,'Timer']],[45,75;130,75;175,95])add_line([sys,'/',['Pulse-Reset',13,'Timer']],[70,40;120,40;175,85])add_line([sys,'/',['Pulse-Reset',13,'Timer']],[70,40;150,40;220,80;245,80])add_line([sys,'/',['Pulse-Reset',13,'Timer']],[270,90;325,90])add_line([sys,'/',['Pulse-Reset',13,'Timer']],[215,95;245,95])%     Finished composite block ['Pulse-Reset',13,'Timer'].set_param([sys,'/',['Pulse-Reset',13,'Timer']],...		'position',[115,59,145,111])%     Subsystem  ['Increase',13,'Detector'].new_system([sys,'/',['Increase',13,'Detector']])set_param([sys,'/',['Increase',13,'Detector']],'Location',[681,484,1074,727])add_block('built-in/Memory',[sys,'/',['Increase',13,'Detector/Memory']])set_param([sys,'/',['Increase',13,'Detector/Memory']],...		'position',[140,155,175,185])add_block('built-in/Switch',[sys,'/',['Increase',13,'Detector/Switch']])set_param([sys,'/',['Increase',13,'Detector/Switch']],...		'Threshold','0.001',...		'position',[285,114,310,146])add_block('built-in/Sum',[sys,'/',['Increase',13,'Detector/Sum6']])set_param([sys,'/',['Increase',13,'Detector/Sum6']],...		'inputs','+-',...		'position',[225,114,240,141])add_block('built-in/Outport',[sys,'/',['Increase',13,'Detector/out_1']])set_param([sys,'/',['Increase',13,'Detector/out_1']],...		'position',[340,120,360,140])add_block('built-in/Inport',[sys,'/',['Increase',13,'Detector/Signal']])set_param([sys,'/',['Increase',13,'Detector/Signal']],...		'position',[50,110,70,130])add_block('built-in/Note',[sys,'/',['Increase',13,'Detector/Increase Detector compares Signal at time t to signal at time  t-1.',13,'When difference is greater than the Switch threshold, it outputs',13,' the difference (positive), 0 otherwise.']])set_param([sys,'/',['Increase',13,'Detector/Increase Detector compares Signal at time t to signal at time  t-1.',13,'When difference is greater than the Switch threshold, it outputs',13,' the difference (positive), 0 otherwise.']],...		'position',[210,20,215,25])add_block('built-in/Constant',[sys,'/',['Increase',13,'Detector/This output means',13,'no detection']])set_param([sys,'/',['Increase',13,'Detector/This output means',13,'no detection']],...		'Value','0',...		'position',[235,160,255,180])add_line([sys,'/',['Increase',13,'Detector']],[260,170;262,170;262,140;280,140])add_line([sys,'/',['Increase',13,'Detector']],[245,130;257,130;257,120;280,120])add_line([sys,'/',['Increase',13,'Detector']],[245,130;280,130])add_line([sys,'/',['Increase',13,'Detector']],[180,170;197,170;197,135;220,135])add_line([sys,'/',['Increase',13,'Detector']],[75,120;97,120;97,170;135,170])add_line([sys,'/',['Increase',13,'Detector']],[315,130;335,130])add_line([sys,'/',['Increase',13,'Detector']],[75,120;220,120])%     Finished composite block ['Increase',13,'Detector'].set_param([sys,'/',['Increase',13,'Detector']],...		'position',[115,147,140,193])%     Subsystem  ['Movement',13,'Detector'].new_system([sys,'/',['Movement',13,'Detector']])set_param([sys,'/',['Movement',13,'Detector']],'Location',[336,521,955,824])add_block('built-in/Outport',[sys,'/',['Movement',13,'Detector/Movement End  Pulse Signal']])set_param([sys,'/',['Movement',13,'Detector/Movement End  Pulse Signal']],...		'Port','2',...		'position',[520,140,540,160])add_block('built-in/Outport',[sys,'/',['Movement',13,'Detector/Movement Count']])set_param([sys,'/',['Movement',13,'Detector/Movement Count']],...		'position',[520,100,540,120])add_block('built-in/Demux',[sys,'/',['Movement',13,'Detector/Demux']])set_param([sys,'/',['Movement',13,'Detector/Demux']],...		'outputs','2',...		'position',[435,105,450,125])add_block('built-in/Inport',[sys,'/',['Movement',13,'Detector/Horizontal',13,'Position']])set_param([sys,'/',['Movement',13,'Detector/Horizontal',13,'Position']],...		'position',[35,30,55,50])add_block('built-in/Mux',[sys,'/',['Movement',13,'Detector/Mux1']])set_param([sys,'/',['Movement',13,'Detector/Mux1']],...		'position',[285,98,315,132])add_block('built-in/Inport',[sys,'/',['Movement',13,'Detector/Veritcal',13,'Position']])set_param([sys,'/',['Movement',13,'Detector/Veritcal',13,'Position']],...		'Port','2',...		'position',[35,85,55,105])add_block('built-in/Mux',[sys,'/',['Movement',13,'Detector/Mux2']])set_param([sys,'/',['Movement',13,'Detector/Mux2']],...		'inputs','2',...		'position',[150,121,180,154])add_block('built-in/Fcn',[sys,'/',['Movement',13,'Detector/Vectorial',13,'Velocity']])set_param([sys,'/',['Movement',13,'Detector/Vectorial',13,'Velocity']],...		'Expr','sqrt(u[1]*u[1]+u[2]*u[2])',...		'position',[200,123,230,157])add_block('built-in/Derivative',[sys,'/',['Movement',13,'Detector/XDerivative']])set_param([sys,'/',['Movement',13,'Detector/XDerivative']],...		'position',[95,110,125,130])add_block('built-in/Derivative',[sys,'/',['Movement',13,'Detector/YDerivative']])set_param([sys,'/',['Movement',13,'Detector/YDerivative']],...		'position',[95,145,125,165])add_block('built-in/Clock',[sys,'/',['Movement',13,'Detector/Clock']])set_param([sys,'/',['Movement',13,'Detector/Clock']],...		'position',[215,185,235,205])add_block('built-in/MATLAB Fcn',[sys,'/',['Movement',13,'Detector/Detector Algorithm']])set_param([sys,'/',['Movement',13,'Detector/Detector Algorithm']],...		'MATLAB Fcn','movedetect(u(1),u(2),u(3),u(4))',...		'Output Width','2',...		'position',[350,73,410,157])add_block('built-in/Note',[sys,'/',['Movement',13,'Detector/Inputs are X,Y, velocity, and time']])set_param([sys,'/',['Movement',13,'Detector/Inputs are X,Y, velocity, and time']],...		'position',[380,175,385,180])add_block('built-in/Note',[sys,'/',['Movement',13,'Detector/Output is the number of movements that have',13,'been completed                     						                                     ']])set_param([sys,'/',['Movement',13,'Detector/Output is the number of movements that have',13,'been completed                     						                                     ']],...		'position',[410,190,415,195])add_block('built-in/Note',[sys,'/',['Movement',13,'Detector/mvmnt[count,1:9]=[amp,dur,pvel,ontime,time,onx,ony,x,y]']])set_param([sys,'/',['Movement',13,'Detector/mvmnt[count,1:9]=[amp,dur,pvel,ontime,time,onx,ony,x,y]']],...		'position',[435,255,440,260])add_block('built-in/Note',[sys,'/',['Movement',13,'Detector/measurements are stored in the global variable "mvmnt"']])set_param([sys,'/',['Movement',13,'Detector/measurements are stored in the global variable "mvmnt"']],...		'position',[430,240,435,245])add_line([sys,'/',['Movement',13,'Detector']],[455,120;475,120;475,150;515,150])add_line([sys,'/',['Movement',13,'Detector']],[455,110;515,110])add_line([sys,'/',['Movement',13,'Detector']],[415,115;430,115])add_line([sys,'/',['Movement',13,'Detector']],[60,40;80,40;90,120])add_line([sys,'/',['Movement',13,'Detector']],[60,95;70,95;70,155;90,155])add_line([sys,'/',['Movement',13,'Detector']],[235,140;240,140;240,120;280,120])add_line([sys,'/',['Movement',13,'Detector']],[60,95;205,95;205,110;280,110])add_line([sys,'/',['Movement',13,'Detector']],[240,195;260,195;260,130;280,130])add_line([sys,'/',['Movement',13,'Detector']],[185,140;195,140])add_line([sys,'/',['Movement',13,'Detector']],[130,120;130,130;145,130])add_line([sys,'/',['Movement',13,'Detector']],[130,155;130,145;145,145])add_line([sys,'/',['Movement',13,'Detector']],[60,40;250,40;250,100;280,100])add_line([sys,'/',['Movement',13,'Detector']],[320,115;345,115])%     Finished composite block ['Movement',13,'Detector'].set_param([sys,'/',['Movement',13,'Detector']],...		'position',[40,148,65,192])%     Subsystem  ['Auto-',13,'Stop'].new_system([sys,'/',['Auto-',13,'Stop']])set_param([sys,'/',['Auto-',13,'Stop']],'Location',[603,361,806,519])add_block('built-in/Inport',[sys,'/',['Auto-',13,'Stop/in_2']])set_param([sys,'/',['Auto-',13,'Stop/in_2']],...		'Port','2',...		'position',[15,80,35,100])add_block('built-in/Inport',[sys,'/',['Auto-',13,'Stop/in_1']])set_param([sys,'/',['Auto-',13,'Stop/in_1']],...		'position',[15,40,35,60])add_block('built-in/Note',[sys,'/',['Auto-',13,'Stop/Stops simulation when the',13,'2 inputs are equal']])set_param([sys,'/',['Auto-',13,'Stop/Stops simulation when the',13,'2 inputs are equal']],...		'position',[95,5,100,10])add_block('built-in/Stop Simulation',[sys,'/',['Auto-',13,'Stop/Stop Simulation']])set_param([sys,'/',['Auto-',13,'Stop/Stop Simulation']],...		'position',[140,50,175,90])add_block('built-in/Relational Operator',[sys,'/',['Auto-',13,'Stop/Equals']])set_param([sys,'/',['Auto-',13,'Stop/Equals']],...		'Operator','==',...		'position',[85,58,115,82])add_line([sys,'/',['Auto-',13,'Stop']],[120,70;135,70])add_line([sys,'/',['Auto-',13,'Stop']],[40,90;60,90;60,75;80,75])add_line([sys,'/',['Auto-',13,'Stop']],[40,50;60,50;60,65;80,65])%     Finished composite block ['Auto-',13,'Stop'].set_param([sys,'/',['Auto-',13,'Stop']],...		'position',[40,59,65,106])add_block('built-in/Note',[sys,'/','My Simulink Library'])set_param([sys,'/','My Simulink Library'],...		'position',[175,10,180,15])%     Subsystem  ['On-Off',13,'Timer'].new_system([sys,'/',['On-Off',13,'Timer']])set_param([sys,'/',['On-Off',13,'Timer']],'Location',[374,266,742,443])add_block('built-in/Clock',[sys,'/',['On-Off',13,'Timer/Clock1']])set_param([sys,'/',['On-Off',13,'Timer/Clock1']],...		'position',[45,30,65,50])add_block('built-in/Memory',[sys,'/',['On-Off',13,'Timer/Memory']])set_param([sys,'/',['On-Off',13,'Timer/Memory']],...		'position',[95,90,130,120])add_block('built-in/Switch',[sys,'/',['On-Off',13,'Timer/Switch']])set_param([sys,'/',['On-Off',13,'Timer/Switch']],...		'Threshold','0.001',...		'position',[180,79,210,111])add_block('built-in/Note',[sys,'/',['On-Off',13,'Timer/Increments time only when the Signal is greater than',13,'Switch threshold.']])set_param([sys,'/',['On-Off',13,'Timer/Increments time only when the Signal is greater than',13,'Switch threshold.']],...		'position',[195,5,200,10])add_block('built-in/Outport',[sys,'/',['On-Off',13,'Timer/time']])set_param([sys,'/',['On-Off',13,'Timer/time']],...		'position',[320,85,340,105])add_block('built-in/Inport',[sys,'/',['On-Off',13,'Timer/on-off']])set_param([sys,'/',['On-Off',13,'Timer/on-off']],...		'position',[45,65,65,85])add_line([sys,'/',['On-Off',13,'Timer']],[215,95;280,95;280,155;75,155;75,105;90,105])add_line([sys,'/',['On-Off',13,'Timer']],[70,40;120,40;175,85])add_line([sys,'/',['On-Off',13,'Timer']],[70,75;130,75;175,95])add_line([sys,'/',['On-Off',13,'Timer']],[135,105;175,105])add_line([sys,'/',['On-Off',13,'Timer']],[215,95;315,95])%     Finished composite block ['On-Off',13,'Timer'].set_param([sys,'/',['On-Off',13,'Timer']],...		'position',[185,59,215,111])%     Subsystem  ['Decrease',13,'Detector'].new_system([sys,'/',['Decrease',13,'Detector']])set_param([sys,'/',['Decrease',13,'Detector']],'Location',[681,484,1074,727])add_block('built-in/Memory',[sys,'/',['Decrease',13,'Detector/Memory']])set_param([sys,'/',['Decrease',13,'Detector/Memory']],...		'position',[140,155,175,185])add_block('built-in/Switch',[sys,'/',['Decrease',13,'Detector/Switch']])set_param([sys,'/',['Decrease',13,'Detector/Switch']],...		'Threshold','1',...		'position',[285,114,310,146])add_block('built-in/Sum',[sys,'/',['Decrease',13,'Detector/Sum6']])set_param([sys,'/',['Decrease',13,'Detector/Sum6']],...		'inputs','-+',...		'position',[225,114,240,141])add_block('built-in/Outport',[sys,'/',['Decrease',13,'Detector/out_1']])set_param([sys,'/',['Decrease',13,'Detector/out_1']],...		'position',[340,120,360,140])add_block('built-in/Inport',[sys,'/',['Decrease',13,'Detector/Signal']])set_param([sys,'/',['Decrease',13,'Detector/Signal']],...		'position',[50,110,70,130])add_block('built-in/Note',[sys,'/',['Decrease',13,'Detector/Decrease Detector compares Signal at time t to signal at time  t-1.',13,'When difference is greater than the Switch threshold, it outputs',13,' the difference (positive), 0 otherwise.']])set_param([sys,'/',['Decrease',13,'Detector/Decrease Detector compares Signal at time t to signal at time  t-1.',13,'When difference is greater than the Switch threshold, it outputs',13,' the difference (positive), 0 otherwise.']],...		'position',[210,20,215,25])add_block('built-in/Constant',[sys,'/',['Decrease',13,'Detector/This output means',13,'no detection']])set_param([sys,'/',['Decrease',13,'Detector/This output means',13,'no detection']],...		'Value','0',...		'position',[235,160,255,180])add_line([sys,'/',['Decrease',13,'Detector']],[260,170;262,170;262,140;280,140])add_line([sys,'/',['Decrease',13,'Detector']],[245,130;257,130;257,120;280,120])add_line([sys,'/',['Decrease',13,'Detector']],[245,130;280,130])add_line([sys,'/',['Decrease',13,'Detector']],[180,170;197,170;197,135;220,135])add_line([sys,'/',['Decrease',13,'Detector']],[75,120;97,120;97,170;135,170])add_line([sys,'/',['Decrease',13,'Detector']],[315,130;335,130])add_line([sys,'/',['Decrease',13,'Detector']],[75,120;220,120])%     Finished composite block ['Decrease',13,'Detector'].set_param([sys,'/',['Decrease',13,'Detector']],...		'orientation',2,...		'position',[195,147,220,193])drawnow% Return any arguments.if (nargin | nargout)	% Must use feval here to access system in memory	if (nargin > 3)		if (flag == 0)			eval(['[ret,x0,str,ts,xts]=',sys,'(t,x,u,flag);'])		else			eval(['ret =', sys,'(t,x,u,flag);'])		end	else		[ret,x0,str,ts,xts] = feval(sys);	endelse	drawnow % Flash up the model and execute load callbackend