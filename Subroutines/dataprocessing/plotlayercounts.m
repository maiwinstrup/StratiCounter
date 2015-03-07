function plotlayercounts(layercounts,data)

%% plotlayercounts(layercounts,data)
% Plot layer counts onto current figure. Uncertain layers are plotted 
% smaller and as dotted lines. Dimensions of line depend on the data range. 
% Copyright (C) 2015  Mai Winstrup

%% Plot layer counts onto figure:
% Bar sizes:
hbar0 = quantile(data,0.05);
hbar1 = quantile(data,0.95);
hbar2 = quantile(data,0.90);  

% Plot layer counts:
if ~isempty(layercounts)
    % Certain layers:
    mask = layercounts(:,3)==0;
    plot(layercounts(mask,1)*[1 1],[hbar0 hbar1],'-r')
    % Uncertain layers:
    if sum(layercounts(:,3))>0
       plot(layercounts(~mask,1)*[1 1],[hbar0 hbar2],'--r')
    end
end