function flow = compute1dFlow(xt_image,blur,dt_filt,dx_filt,sig)
% flow = compute1dFlow(xt_image,blur,dt_filt,dx_filt,sig)
%
% Compute flow of an XT image.  BLUR is a 2d filter kernel (optional).
% DT is a 2D temporal derivative kernel.

if (exist('blur') ~= 1)
	blur = [0.25 0.5 0.25]'*[0.25 0.5 0.25];
end

if (exist('dt_filt') ~= 1)
	dt_filt  = [-0.459 0.0 0.459]' * [0.233 0.534 0.233];
end

if (exist('dx_filt') ~= 1)
	dx_filt  = [0.233 0.534 0.233]' * [-0.459 0.0 0.459];
end

if (exist('sig') ~= 1)
	sig  = 1e-4;
end

dx_im = corrDn(xt_image,dx_filt,'dont-compute');
dt_im = corrDn(xt_image, dt_filt,'dont-compute');

num = corrDn(dt_im.*dx_im, blur, 'dont-compute');
den = corrDn(dx_im.*dx_im, blur, 'dont-compute');

den = den + sig;

flow = - num ./ den;
