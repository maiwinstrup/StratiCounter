function [manualcounts, meanLambda, newinterval] = loadlayercounts(Model,interval)

%% [manualcounts, meanLambda, newinterval] = loadlayercounts(Model,interval)
% Load manually counted layers for the depth interval in consideration. 
% If manual layer counts do not exist for the entire interval, they are 
% provided for "newinterval" only. 

% Mai Winstrup
% 2014-06-17 21:11: Updates to input and output of function
% 2014-08-21 20:08: Major restructuring, a data file containting the manual 
%                   counts are now assumed to exist beforehand.
% 2014-08-22 16:30: Changes to what is "outside interval", ensure new
%                   interval fits to our depthscale
% 2014-08-23 13:20: Ensure layer counts are sorted
% 2014-10-08 09:58: Layer counts are located halfway between data points.

%% Import manual layer counts:
manualcounts = importdata(Model.pathManualCounts);
if isstruct(manualcounts); manualcounts = manualcounts.data; end
% Ensure that they are sorted according to increasing depth: 
manualcounts = sortrows(manualcounts,1);

%% Remove layercounts from outside interval, and ensure layer positions to 
% be located halfway between data points:
% Location of first and last data point:
dstart = ceil((interval(1)-Model.dx_center)/Model.dx)*Model.dx+Model.dx_center; 
dend = ceil((interval(2)-Model.dx_center)/Model.dx)*Model.dx+Model.dx_center;
% Possible layer boundary locations:
depth_px = dstart-Model.dx/2:Model.dx:dend+Model.dx/2;

% Truncated locations:
pos = interp1(depth_px,1:length(depth_px),manualcounts(:,1),'nearest',nan);
% Layers outside interval are nan, these are removed: 
mask = isfinite(pos);
manualcounts = manualcounts(mask,:);

% If no manual layer counts exist in interval:
% Return to main program with error message.
if isempty(manualcounts)
    disp(['Manual counts are not provided for interval ' num2str(interval(1)) ...
        '-' num2str(interval(2)) 'm. Please correct!'])
    return
end

% Replacing with truncated layer boundary locations in manualcounts:
manualcounts(:,1) = depth_px(pos(mask));

%% Mean layer thickness in interval:
meanLambda = mean(diff(manualcounts(:,1)));

%% If layer counts do exist for part of the interval, but not all:
% The interval boundaries are changed to reflect this (missing layer counts 
% in a small area close to the boundaries are allowed).
newinterval = interval;
if manualcounts(1,1)>interval(1)+5*meanLambda
    newinterval(1) = manualcounts(1,1);
end
if manualcounts(end,1)<interval(2)-5*meanLambda;
    newinterval(2) = manualcounts(end,1);
end