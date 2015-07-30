function [stack, dlayer] = stacklayers(depth,data,layercounts,unc_layercounts,dt)

%% [stack, dlayer] = stacklayers(depth,data,layercounts,unc_layercounts,dt):
% Stacking data as annual layers. All layers are linearly interpolated to 
% have the same length, with each data point corresponding to timestep dt. 
% dlayer contains the starting depth of the stacked layers.
% If uncertain layers exist, a composite stack is obtained from two stacks 
% in which every second uncertain layer is removed. Certain layers will
% hence occur twice in the final stack, while uncertain layers occur only 
% once. 
% Copyright (C) 2015  Mai Winstrup

%% Removing data from outside the first and last certain annual layer mark:
if ~isempty(unc_layercounts)
    % Remove uncertain layers from start and end:
    istart = find(unc_layercounts==0,1,'first');
    iend = find(unc_layercounts==0,1,'last');
    layercounts = layercounts(istart:iend);
end
% Removing data from outside first and last layer boundary:
mask = depth>=layercounts(1,1)-10^-10 & depth<=layercounts(end,1)+10^-10;
depth = depth(mask);
data = data(mask);

%% Stack the data record in years: 
% If no uncertain layers exist: 
if sum(unc_layercounts)==0
    % Stack all layers:
    stack = makestack(depth,data,layercounts,dt);
    dlayer = layercounts(1:end-1);
    
else
    % Make a composite stack consisting of the two stacks obtained by 
    % removing every second uncertain layer boundary. 
    
    % Uncertain layer boundaries:
    iunc = find(unc_layercounts==1);
    
    % Stack1: Removing every second of these:
    mask1 = true(size(layercounts));
    mask1(iunc(1:2:end))=false; % These are removed
    % Make stack:
    stack1 = makestack(depth,data,layercounts(mask1),dt);
    dlayer1 = layercounts(mask1);
    dlayer1 = dlayer1(1:end-1);
    
    % Stack2: Same but keeping the other half of uncertain boundaries: 
    mask2 = true(size(layercounts));
    mask2(iunc(2:2:end))=false; % These are removed
    stack2 = makestack(depth,data,layercounts(mask2),dt);
    dlayer2 = layercounts(mask2);
    dlayer2 = dlayer2(1:end-1);

    % Composite stack:
    stack = [stack1; stack2];
    dlayer = [dlayer1; dlayer2]; % Start of layers in stack
end
end

%% Embedded interpolation and stacking function: 
function stack = makestack(d,y,dlayer,dt)

%% stack = makestack(d,y,dlayer,dt)
% Stack the data (d,y) according to the layer boundary positions (dlayer).
% These boundaries are taken as the location midway between two data points
% belonging to different layers. 
% Mai Winstrup

%% Remove data from outside interval:
mask = d>=dlayer(1)&d<=dlayer(end);
d = d(mask);
y = y(mask);

% Add a NaN data point before and after (required for interpolation):
d = [d(1)-(d(2)-d(1)); d(:); d(end)+(d(end)-d(end-1))];
y = [NaN;y(:);NaN];

%% Up/downsample data to given dt:
% Construct preliminary timescale for data: 
% Assuming all layers equals 1 year, and using linear interpolation 
% inbetween. 
% Age for each data point:
tpx = interp1(dlayer,0:1:(length(dlayer)-1),d,'linear','extrap');

% Up/downsample data to given dt:
% Do not consider layer fractions at edges.
tnew = 1/2*dt:dt:floor(tpx(end)); % Times are taken half way
ynew = downsampling(tpx,y,tnew);

%% Stack data:
stack = reshape(ynew,1/dt,length(ynew)*dt)';
end