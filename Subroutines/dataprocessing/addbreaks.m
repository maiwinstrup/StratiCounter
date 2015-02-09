function [depth_new, data_new, breaks] = addbreaks(depth,data,dx) 

%% [depth_new, data_new, breaks] = addbreaks(depth,data,dx);
% Find locations where the difference between measurements is greater than
% dx. These are areas of breaks. For these sections, depths are added to
% file, along with a corresponding NaN data value. 
% Copyright (C) 2015  Mai Winstrup

%% Find areas where differences between measurements is larger than dx.
diffdepth = abs(diff(depth));
index = find(diffdepth>dx);

%% Add extra depth entries onto the end:
depths_nan = depth(index) + dx;
breaks = depths_nan; % Starting depth of breaks
data_nan = NaN(size(depths_nan));
depth_new = [depth; depths_nan];
data_new = [data; data_nan];

%% Sorting them: 
[depth_new, index] = sort(depth_new);
data_new = data_new(index);