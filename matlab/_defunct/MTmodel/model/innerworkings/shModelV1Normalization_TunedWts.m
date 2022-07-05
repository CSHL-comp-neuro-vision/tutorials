function [wts, kori, ksf, ktf] = shtunenrmwts(newDirs, dirs)

% for i = 1:size(newDirs, 1)
%     for j = 1:size(dirs, 1)
%         
%         newdir = newDirs(i,:);
%         dir = dirs(j,:);
% 
%         pref1 = v12sin(newdir);
%         newori = pref1(1);
%         newsf = pref1(2);
%         newtf = pref1(3);
%         
%         pref2 = v12sin(newdir);
%         ori = pref2(1);
%         sf = pref2(2);
%         tf = pref2(3);
%         
%    
%         orisig = 3/5.*pi./2;
%         sfsig = log2(8);
%         tfsig = log2(8);
%         
%         k = pi - ori;
%         ori = ori + k;
%         newori = newori + k;
%         newori = mod(newori, 2*pi);
%         
%         xori = mod(abs(ori - newori), pi);
%           
%         xsf = sf./newsf;
%         xtf = tf./newtf;
%         
%         kori = 1./sqrt(2.*pi.*orisig.^2) .* exp(-(xori).^2./(orisig.^2));
%         ksf = 1./sqrt(2.*pi.*sfsig.^2) .* exp(-(log2(xsf).^2)./(sfsig.^2));
%         ktf = 1./sqrt(2.*pi.*tfsig.^2) .* exp(-(log2(xtf).^2)./(tfsig.^2));
%         
%         wts(i,j) = kori.*ksf.*ktf;
%     end
% end
% 
% 
% n1 = 1./sqrt(2.*pi.*orisig.^2);
% n2 = 1./sqrt(2.*pi.*sfsig.^2);
% n3 = 1./sqrt(2.*pi.*tfsig.^2);
% wts = wts./(n1.*n2.*n3).*.36./.084;

wts = ones(size(newDirs, 1), size(dirs, 1));
