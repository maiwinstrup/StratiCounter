function [dstart_fig,dend_fig,depth_fig,data_fig] = ...
    depthrangefigure(depth,data,layercounts)

%% [dstart_fig,dend_fig,depth_fig,data_fig] = 
%       depthrangefigure(depth,data,layercounts)
% Select data section to display on figure. Data section is selected to 
% contain approximately 25 layers, and NaNs in data series should cover 
% less than 50% of section. 
% Copyright (C) 2015  Mai Winstrup

%% Depth section to be shown in plot: 
% Select section starting from the beginning containing approximately 25
% layers:
dstart_fig = depth(1);
if ~isempty(layercounts)
    L = mean(diff(layercounts(:,1)))*25; % Approximately 25 layers   
else
    L = 1; % If layer counts do not exist: Display a 1 m section. 
end
% Corresponding end of data section:
dend_fig = min(dstart_fig+L,depth(end));
mask = depth>=dstart_fig & depth<=dend_fig;
% Depthscale for figure:
depth_fig = depth(mask);
data_fig = data(mask);

% Data should cover at least half of this section:
% If not, we move to the end of current section and try to start another 
% section here. 
while sum(isnan(data_fig))>0.5*length(depth_fig)
    dstart_fig = dend_fig;
    dend_fig = min(dstart_fig+L,depth(end));
    mask = depth>=dstart_fig & depth<=dend_fig;
    depth_fig = depth(mask);
    data_fig = data(mask);
    
    % In case this is never possible: 
    % Display data starting from the beginning of data series.   
    if dend_fig==depth(end)
        dstart_fig = depth(1);
        dend_fig = depth(1)+L;
        break; 
    end
end