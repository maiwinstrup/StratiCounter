function ynew = smoothdata(x,y,L)

%% ynew = smoothdata(x,y,L)
% Smooth data over a section L, measured in same units as x. 
% Smoothing is performed using a moving average.
% Copyright (C) 2015  Mai Winstrup

%% Upsample data to equidistant data values if necessary:
% If equidistant data, the upsampling and downsampling here will have no 
% effect. 
if mean(diff(x))>=min(diff(x))
    dxmin = min(diff(x));
    x_up = x;
    y_up = y;
else
    dxmin = nanmin(diff(x))/5; 
    x_up = x(1):dxmin:x(end);
    y_up = interp1(x,y,x_up);
end

%% Smooth upsampled data:
yup_smooth = smooth(y_up,round(L/dxmin));

%% Downsample to original resolution:
ynew = interp1q(x_up(:),yup_smooth(:),x(:));