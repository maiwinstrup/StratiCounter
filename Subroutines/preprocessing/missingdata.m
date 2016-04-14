function [depth_out, data_out] = missingdata(depth_in,data_in) 

%% [depth_out, data_out] = missingdata(depth_in,data_in)
% Find missing sections in data series. For these sections, a depth value 
% corresponding to start and end of the break is added to the file, along 
% with two corresponding NaN data values. 
% Copyright (C) 2015  Mai Winstrup

%% Set threshold value: 
% Where spacing are significantly larger than the 99% percentile of the 
% distribution. 
dx_in = diff(depth_in);
% Running threshold value:
threshold = ones(length(depth_in),1); 
for i = 101:200:length(depth_in)-99;
    % Values are taken for every 200 data point, with area spanning 400 
    % data points total
    istart = max(1,i-200);
    iend = min(length(depth_in)-1,i+200);
    threshold(i-100:i+99)=1.5*quantile(dx_in(istart:iend),0.99);
end
% Extend with similar values:
threshold(i+1:length(depth_in))=threshold(i);

%% Find sections where differences between measurements is larger than 
% threshold:
diffdepth = abs(diff(depth_in));
index = find(diffdepth>threshold(1:end-1));

%% Add extra depth entries onto the end:
depths_nan1 = depth_in(index) + median(dx_in);
depths_nan2 = depth_in(index+1) - median(dx_in);
depths_nan = [depths_nan1; depths_nan2];
startofbreaks = depths_nan1 - median(dx_in); % Starting depth of breaks
endofbreaks = depth_in(index+1); % Ending depth of breaks

data_nan = NaN(length(depths_nan),size(data_in,2),size(data_in,3));
depth_out = [depth_in; depths_nan];
data_out = [data_in; data_nan];

%% Sorting them: 
[depth_out, index] = sort(depth_out);
data_out = data_out(index,:,:);

%% Display information:
if ~isempty(startofbreaks)
    disp('Break(s) are added to datafile starting at the following depth(s): ');
    disp(num2str(startofbreaks))
    disp('and ending at depth(s): ')
    disp(num2str(endofbreaks))
end