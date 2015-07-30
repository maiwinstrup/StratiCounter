function [depth_new, data_new, startofbreaks, endofbreaks] = addbreaks(depth,data,dx) 

%% [depth_new, data_new, startofbreaks, endofbreaks] = addbreaks(depth,data,dx)
% Find locations where the difference between measurements is greater than
% dx. These are areas of breaks. For these sections, a depth value 
% corresponding to start of the break is added to the file, along with a 
% corresponding NaN data value. 
% Copyright (C) 2015  Mai Winstrup

%% Find areas where differences between measurements is larger than dx.
diffdepth = abs(diff(depth));
index = find(diffdepth>dx);

%% Add extra depth entries onto the end:
depths_nan = depth(index) + dx;
startofbreaks = depths_nan; % Starting depth of breaks
endofbreaks = depth(index+1); % Ending depth of breaks

data_nan = NaN(length(depths_nan),size(data,2),size(data,3));
depth_new = [depth; depths_nan];
data_new = [data; data_nan];

%% Sorting them: 
[depth_new, index] = sort(depth_new);
data_new = data_new(index,:,:);