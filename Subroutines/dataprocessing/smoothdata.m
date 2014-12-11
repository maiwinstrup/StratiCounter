function ynew = smoothdata(x,y,L)

%% ynew = smoothdata(x,y,L)
% Smooth data over a section L, measured in same units as x. 
% Smoothing is performed using a moving average.
% Mai Winstrup, 2014
% 2014-07-16 13:58

%% Upsample data to equidistant data values:
% If equidistant data, the upsampling and downsampling here will have no 
% effect. 
dxmin = nanmin(diff(x))/5; 
x_up = x(1):dxmin:x(end);
y_up = interp1(x,y,x_up);

%% Smooth upsampled data:
y_up_smooth = smooth(y_up,L/dxmin);

%% Downsample to original resolution:
ynew = interp1q(x_up(:),y_up_smooth(:),x(:));