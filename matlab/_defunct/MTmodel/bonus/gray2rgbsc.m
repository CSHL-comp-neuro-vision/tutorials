% MRGB = gray2rgb(M) - put a 2D matrix into RGB format for image display
% If M is a 2D matrix, MRGB is an RGB image that will look
% grayscale and can be displayed using image. The minimum of M is
% first set to zero, and M is pointwise divided so that its maximum
% value is 1.

function mrgb = gray2RGB(m)

if size(size(m), 2) == 2
    m = m - min(min(m));
    m = m./max(max(m));
    mrgb(:,:,1) = m;
    mrgb(:,:,2) = 0;
    mrgb(:,:,3) = 0;
    mrgb = ntsc2rgb(mrgb);
    
else
    fprintf('gray2rgb works only on 2D matrices. So sorry.');
    mrgb = [];
end