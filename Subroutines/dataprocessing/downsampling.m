function [depth_new, data_new] = ...
    downsampling(depth,data,dx_new,dx_center,dstart,dend,method,plotlevel)

%% [depth_new, data_new] = downsampling(depth,data,dx_new,dx_midpoint,method,plotlevel)
% Downsampling using a step function:
% Can be done in two ways: 
% Method 1: Upsampling, smoothing with boxcar function, and subsequently
% downsampling. Values calculated only for areas fully covered (the 
% default, is faster).
% Method 2: Data values are calculated independently for areas with more
% than 50% coverage in original data (slow). 
% Both: Data are allowed not to be on equidistant depth scale.

% Copyright (C) 2015  Mai Winstrup
% 2014-04-21 17:08: Changed calculation for constantsampling_test and
%                   samplerate_test
% 2014-05-18 19:41: Minor adjustments
% 2014-07-16 15:19: showplots->plotlevel, 
%                   using dx_midpoint to determine new depth scale
% 2014-07-24 10:57: dstart found using floor instead of round
% 2014-08-06 17:15: Correction to find of new depthscale: using round 
%                   instead of floor
% 2014-08-21 22:33: -do-

%% Setting default values:

if nargin < 4; dx_center = 0; end
if nargin < 5; dstart = depth(1); end
if nargin < 6; dend = depth(end); end
if nargin < 7; method = 'method1'; end
if nargin < 8; plotlevel = 0; end

%% New depthscale:
% Starting point for resulting depthscale:
dstart = ceil((dstart-dx_center)/dx_new)*dx_new+dx_center;
dend = ceil((dend-dx_center)/dx_new)*dx_new+dx_center;
depth_new(:,1) = dstart:dx_new:dend;
% Round to increase accuracy:
depth_new(:,1) = round((depth_new(:,1)-dx_center)/dx_new)*dx_new+dx_center;

%% Case 1: Data already has the right sampling distance, and correct 
% starting point. 
% In thise case, we simply interpolate to make sure the data points are
% placed appropriately. Further, eventual missing sections (sections without 
% corresponding depth entry) are added, and given data values of NaNs.

% Is data already equidistant?
constantsampling_test = abs(min(diff(depth))-max(diff(depth)))/mean(diff(depth));
% And with correct sample rate?
samplerate_test = abs(mean(diff(depth))-dx_new)/mean(diff(depth));
% Check for correct starting point:
diffstartpoint = abs(dstart-depth(1))/dx_new;
% Measured as fraction of the distance between data points in resulting 
% timescale. 

if constantsampling_test<0.5 && samplerate_test<10^-3 && diffstartpoint<0.1; 
    % Check for sections of nan, and add in NaNs if found:
    % These are selected as sections with depth differences larger than
    % 1.5*dx_new.
    [depth, data] = addbreaks(depth,data,1.5*dx_new);
    
    % Interpolate to new depth scale:
    data_new = interp1q(depth,data,depth_new);
    return
end

%% Case 2: Downsampling to equidistant depth scale:
% If we have a constant sample rate, we check for sections of breaks in 
% data file: 
if constantsampling_test < 0.5
    % Check for breaks, sections which do not have a corresponding depth 
    % entry, and add in NaNs if found.
    % Areas of breaks are selected as sections with depth differences 
    % larger than 1.5*dx_old.
    dx_old = mean(diff(depth));
    [depth, data] = addbreaks(depth,data,1.5*dx_old);
end

switch method
    case 'method1'    
        %% Upsampling the data to e.g. 100 points per new value of depth:
        dx_upsample = dx_new/100;
        depth_upsample(:,1) = floor(depth(1)/dx_upsample)*dx_upsample:dx_upsample:ceil(depth(end)/dx_upsample)*dx_upsample;
        data_upsample = interp1q(depth,data,depth_upsample);
        
        % Filtering using a box-car function:
        b = 1/101*ones(1,101);
        data_filt = filter(b,1,data_upsample);
        data_filt = [data_filt(51:end); nan(50,1)];
         
        % Downsampling the filtered data using interpolation:        
        data_new = interp1q(depth_upsample,data_filt,depth_new);
        
    case 'method2'
        % Number of new data values:
        N = length(depth_new);

        %% Averaging data within these sections:
        % Boundaries corresponding to current dataervations:
        depth = depth(:);
        data = data(:);
        dbounds = [depth(1); depth(1:end-1)+diff(depth)/2; depth(end)];
        mask = isfinite(data);
        data_new = nan(N,1);

        for i = 1:N
            % Interval boundaries:
            dstart = depth_new(i)-0.5*dx_new;
            dend = depth_new(i)+0.5*dx_new;
        
            % Weights of individual data points:
            t1 = max(dbounds(1:end-1),dstart);
            t2 = min(dend,dbounds(2:end));
            w = max(t2-t1,0); % Coverage of selected section [m]
    
            % Corresponding mean value, provided more than half of section
            % is covered by the original data:
            coverage = sum(mask(w>0).*w(w>0)); 
            if coverage > 0.5*dx_new
                w = w(mask)/sum(w(mask));
                data_new(i)=sum(w.*data(mask));
            end
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