function [dstart_fig,dend_fig, depth_fig,data_fig] = selectdepthintervalforfigure(depth,data,layercounts)
%% Depth section to be shown in plot: 
dstart_fig = depth(1);
if ~isempty(layercounts)
    L = mean(diff(layercounts(:,1)))*25; % Approximately 25 layers   
else
    L = 1;
end
dend_fig = min(dstart_fig+L,depth(end));
mask = depth>=dstart_fig & depth<=dend_fig;
% Depthscale for figure:
depth_fig = depth(mask);
data_fig = data(mask);

% Data should cover at least half of section:
while sum(isnan(data_fig))>0.5*length(depth_fig)
    dstart_fig = dend_fig;
    dend_fig = min(dstart_fig+L,depth(end));
    mask = depth>=dstart_fig & depth<=dend_fig;
    depth_fig = depth(mask);
    data_fig = data(mask);
    
    % If this never is possible: Use data starting from the beginning. 
    if dend_fig==depth(end)
        dstart_fig = depth(1);
        dend_fig = depth(1)+L;
        break; 
    end
end