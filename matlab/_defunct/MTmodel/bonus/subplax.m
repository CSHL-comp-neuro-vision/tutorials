function h = subplax(Nax, ll, wd, flag)
% h = subplax(NumberOfAxes_Y_X, LowerLeft, Ht_Wdth, Flag);
%
% Creates axes for subplots (more powerful and flexible than subplot)
%
% Inputes:
% NumberOfAxes_Y_X = array with [# y plots, # x plots] in grid
% LowerLeft = [y,x] position of lower left corner
% Ht_Wdth = [# of verticle grid blocks, # of horiz grid blocs];
% Flag - determines how much whitespace between axis box and grid
%        0 - matches spacing between matlab "subplot" function
%        1 - exact grid spacing
%        .45 = .45 % spacing evenly distributed around figure (max .99)
%       [lf btm rt top] = % spacing on each side specified
%
%  Numbering of grid looks like 
%    Y =  3    X = 1 2 3 ...
%         2
%         1
%  Note: indices not constrained to be integers.


if nargin < 4
    flag = 0;
end


nrows = Nax(1);
ncols = Nax(2);

if (flag(1) == 0) & (length(flag) == 1)      % Use matlab default spacing
    row = ll(1)-1;
    col = ll(2)-1;    
    Pct_L = 2*0.09;  % percent offsets on each side
    Pct_R = 2*0.045;
    Pct_B = Pct_L;
    Pct_T = Pct_R;
    if Nax(1) > 2
        Pct_T = 0.9*Pct_T;     Pct_B = 0.9*Pct_B;
    end
    if Nax(2) > 2
        Pct_L = 0.9*Pct_L;     Pct_R = 0.9*Pct_R;
    end

    dflt = get(gcf,'DefaultAxesPosition');  % default figure position
    dcol = dflt(3)*(Pct_L+Pct_R)/(ncols-Pct_L-Pct_R); % spacing from lower left corner
    drow = dflt(4)*(Pct_B+Pct_T)/(nrows-Pct_B-Pct_T);
    
    wid = dflt(3) + dcol;
    hht = dflt(4) + drow;
    width = wid/ncols-dcol;
    height = hht/nrows-drow;
    
    xpos = dflt(1)+col*wid/ncols;
    ypos = dflt(2)+row*hht/nrows;
    xpos2 = xpos+(wd(2)-1)*wid/ncols + width;
    ypos2 = ypos+(wd(1)-1)*hht/nrows + height;
    pos = [xpos ypos xpos2-xpos ypos2-ypos];
    
else
    row = 1./Nax(1)*(ll(1)-1);
    col = 1./Nax(2)*(ll(2)-1);
    dy = 1./Nax(1)*wd(1);
    dx = 1./Nax(2)*wd(2);
    
    if flag(1) == 1 % Use subdivision of figure
        pos = [col row dx dy];
    elseif (length(flag) == 1)
        flag = flag/2;
        pos = [col+dx*flag row+dy*flag dx*(1-2*flag) dy*(1-2*flag)];
    elseif (length(flag) == 4) % offsets:  [lft bttm rt top]
        pos = [col+dx*flag(1) row+dy*flag(2) dx*(1-flag(1)-flag(3)) dy*(1-flag(2)-flag(4))];
    else
        error('Flag improperly specified -- subplax.m');
    end
end

h = subplot('position', pos);

