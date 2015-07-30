function [manualcounts, meanLambda, newinterval] = loadlayercounts(Model,interval)

%% [manualcounts, meanLambda, newinterval] = loadlayercounts(Model,interval)
% Load manually counted layers for the depth interval in consideration. 
% If manual layer counts do not exist for the entire interval, they are 
% provided for "newinterval" only. 
% Copyright (C) 2015  Mai Winstrup

%% Import manual layer counts:
manualcounts = importdata(['./Manualcounts/' Model.nameManualCounts]);
if isstruct(manualcounts); manualcounts = manualcounts.data; end

% Ensure that these counts are sorted according to increasing depth: 
manualcounts = sortrows(manualcounts,1);

%% Remove layercounts from outside interval, and ensure correct placement. 
% Layer positions must always be located halfway between data points. 
d = makedepthscale(interval(1),interval(2),Model.dx,Model.dx_offset);
% Possible layer boundary locations:
dlayer_px = [d(:)-Model.dx/2; d(end)+Model.dx/2];
% Truncated locations:
pos = interp1(dlayer_px,1:length(dlayer_px),manualcounts(:,1),'nearest',nan);

% Layers outside interval are nan, these are removed: 
mask = isfinite(pos);
manualcounts = manualcounts(mask,:);

% If no manual layer counts exist in interval:
% Return to main program with error message.
if isempty(manualcounts) 
    warning(['Manual counts are not provided for interval ' num2str(interval(1)) ...
        '-' num2str(interval(2)) 'm'])
    meanLambda = nan; 
    newinterval = interval;
    return
elseif size(manualcounts,1)<=5
    % Warning if very few layer counts exist:
    warning(['Very few manual layer counts ('  num2str(size(manualcounts,1)) ...
        ') exist for interval ' num2str(interval(1)) ...
        '-' num2str(interval(2)) 'm'])
end

% Replacing with truncated layer boundary locations in manualcounts:
manualcounts(:,1) = dlayer_px(pos(mask));

%% Mean layer thickness in interval:
meanLambda = mean(diff(manualcounts(:,1)));

%% If layer counts do exist for part of the interval, but not all:
% The interval boundaries are changed to reflect this (missing layer counts 
% in a small area close to the boundaries is allowed).
newinterval = interval;
if manualcounts(1,1)>interval(1)+5*meanLambda
    newinterval(1) = manualcounts(1,1);
end
if manualcounts(end,1)<interval(2)-5*meanLambda;
    newinterval(2) = manualcounts(end,1);
end