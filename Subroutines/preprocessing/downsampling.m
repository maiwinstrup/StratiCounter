function data_new = downsampling(depth,data,depth_new,plotlevel)

%% data_new = downsampling(depth,data,depth_new,plotlevel)
% Downsampling using a step function: Data is upsampled, smoothing with 
% boxcar function, and subsequently downsampled. 
% Values are subsequently calculated for areas containing nan in data, such
% that a value for these is provided if non-nan data covers more than 50% 
% of an interval. 
% Data on non-equidistant depth scale is allowed.
% Copyright (C) 2015  Mai Winstrup

%% Setting default values:
if nargin < 4; plotlevel = 0; end
 
%% Case 1: Data already has the right sampling distance, and correct 
% starting point. In thise case, we simply interpolate to make sure the 
% data points are placed appropriately. 

% Is data already equidistant?
dx_old = diff(depth);
constantsampling_test = abs(median(dx_old)-min(dx_old))/mean(dx_old);
% And with correct sample rate?
dx_new = mean(diff(depth_new));
samplerate_test = abs(median(dx_old)-dx_new)/dx_new;
% Check for correct values of dx_offset, i.e. must have same first value:  
diffstartpoint = abs(depth_new(1)-depth(1))/dx_new;
% Measured as fraction of the distance between data points in resulting 
% timescale. 

if constantsampling_test<0.01 && samplerate_test<0.01 && diffstartpoint<10^-3; 
    % Interpolate to new depth scale:
    data_new = interp1q(depth,data,depth_new);
    return
end

%% Downsampling to equidistant depth scale:
% Upsampling the data to 100 points per new value of depth:
dx_upsample = dx_new/100;
depth_upsample(:,1) = floor(depth(1)/dx_upsample)*dx_upsample:...
    dx_upsample:ceil(depth(end)/dx_upsample)*dx_upsample;
data_upsample = interp1q(depth,data,depth_upsample);
        
% Average by filtering using a box-car function:
b = 1/101*ones(1,101);
data_filt = filter(b,1,data_upsample);
data_filt = [data_filt(51:end); nan(50,1)];
        
% Downsampling the filtered data using interpolation:        
data_new = interp1q(depth_upsample,data_filt,depth_new(:));
        
% Deal with areas of nan:
nandata_index = find(isnan(data_new));
for i = 1:length(nandata_index)
    % Interval boundaries:
    dstart = depth_new(nandata_index(i))-0.5*dx_new;
    dend = depth_new(nandata_index(i))+0.5*dx_new;

    % Upsampled data section for this interval:
    mask = depth_upsample >= dstart & depth_upsample <= dend;
    datasection_upsample = data_upsample(mask);
            
    % Does data exist for more than 50% of the interval?
    ndatapoints = sum(isfinite(datasection_upsample));
    if ndatapoints>=0.5*length(datasection_upsample)
        % Provide mean value using the 50% coverage of data:
        data_new(nandata_index(i))=nanmean(datasection_upsample);
    end
end

%% Plotting:
if plotlevel>0
    % Making a step curve:
    depth_step = [depth_new(1:end)'-dx_new/2;depth_new(1:end)'+dx_new/2];
    data_step = [data_new';data_new'];
         
    figure;
    plot(depth,data,'.-b')
    hold on
    plot(depth_new,data_new,'xg')
    plot(depth_step(:),data_step(:), '-g', 'linewidth',1.5)
    xlim([depth(1), min(depth(1)+1,depth(end))])
end