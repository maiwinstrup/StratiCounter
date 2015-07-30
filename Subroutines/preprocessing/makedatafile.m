function [data_out,derivnoise,hfigpreproc,hfigderiv] = ...
    makedatafile(data_in,depth_in,preprocsteps,derivatives,depth_out,...
    plotlevel,species,layercounts)

%% [data_out,derivnoise,hfigpreproc,hfigderiv] = makedatafile(data_in,...
%      depth_in,preprocsteps,derivatives,depth_out,plotlevel,species,layercounts)
%
% This function takes in a data series, which may contain multiple species, 
% and preprocesses these according to the specifications in preprocsteps. 
% Subsequently, the data series may be downsampled to the depthscale 
% "depth_out", and derivatives are calculated. "Derivnoise" is a vector 
% providing the theoretical amount of white noise in the derivative data 
% series relative to the original data profile. Layercounts are only used 
% for plotting. Preprocessed data and derivatives are plotted if 
% plotlevel>0. Figure handles are provided as output. 

% Copyright (C) 2015  Mai Winstrup
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.

%% Set default: 
if nargin<8; layercounts=[]; end % No plotting of layer counts
if nargin<7; species=[]; end; % No titles on plot
if nargin<6; plotlevel=0; end % No plotting
if nargin<5; depth_out = depth_in; end % No change in depth scale
    
% If not provided, the species name array is empty:
nSpecies = size(data_in,3);
if isempty(species)
    for j = 1:nSpecies; species{j}=''; end
end

%% Check for unrecorded breaks in data:
% A (likely) break in data can be observed as much larger sections without 
% a data point than "usual". However, this check is only performed if the 
% data points are regularily spaced. If so, we select as breaks those 
% no-data sections which are wider than 1.8x the median value. One depth 
% value (and a data value of nan) corresponding to this break is added to 
% the data records.
    
% Test for regularly spaced data: 
dx_old = diff(depth_in);
regularity_test = abs(median(dx_old)-min(dx_old))/mean(dx_old);
if regularity_test<0.01
    [depth_in, data_in, startofbreaks, endofbreaks] = ...
        addbreaks(depth_in,data_in,1.8*median(dx_old));
    if ~isempty(startofbreaks)
        disp(['Break(s) are added to datafile starting at the following depth(s): ' ...
            num2str(startofbreaks)])
        disp(['and ending at depth(s): ' num2str(endofbreaks)])
    end
end

%% Preprocess data:
data1 = nan(size(data_in));
hfigpreproc = gobjects(nSpecies,1);

for j = 1:nSpecies
    % Check if need to preprocess data: 
    if sum(isfinite(data_in(:,1,j)))==0 || isempty(preprocsteps{j})
        % No need for preprocessing, since:
        % a) it's all NaN anyway, or
        % b) no preprocessing step is specified.
        data1(:,1,j) = data_in(:,1,j);
    else
        % Preprocess data:
        [data1(:,1,j),hfigpreproc(j)] = ...
            preprocessdata(data_in(:,1,j),depth_in,preprocsteps{j},...
            plotlevel,species{j},layercounts);
    end
end

%% Interpolate to depth scale:
% Check for no changes in depthscale
if isequal(depth_out,depth_in); 
    nochangeindepth = true;
else nochangeindepth = false;
end

if nochangeindepth
    data_out(:,1,:)=data1(:,1,:);
    % Initialize rest of array:
    data_out(:,2:1+derivatives.nDeriv,:)=nan;

else
    % Interpolate to new depth scale:
    % Species #1:
    data_out(:,1,1) = downsampling(depth_in,data1(:,1,1),depth_out);    
    % Initialize data array for remaining species (and derivatives):
    L = length(depth_out);
    data_out = cat(3,[data_out(:,1,1),nan(L,derivatives.nDeriv)], ...
        nan(L,1+derivatives.nDeriv,nSpecies-1));
    
    % Remaining species:
    for j = 2:nSpecies
        data_out(:,1,j) = downsampling(depth_in,data1(:,1,j),depth_out);
    end
end

%% Calculate derivatives:
% The derivatives are calculated per pixel, and thus their calculation 
% do not require an equidistant timescale.
hfigderiv = gobjects(nSpecies,1);
derivnoise = [];
for j = 1:nSpecies
    if sum(isfinite(data_out(:,1,j)))==0 
        % Unnecessary if: 
        % a) all data is nan
        %    In this case, all derivatives are nan too (as initialized)        

    elseif nochangeindepth && isempty(preprocsteps{j})
        % b) already calculated using this depthscale, i.e. no processing 
        %    or depth interpolation was performed on data.
        data_out(:,2:1+derivatives.nDeriv,j) = data_in(:,2:1+derivatives.nDeriv,j);
    
    else
        % Otherwise, calculate derivatives:
        [slope,derivnoise,hfigderiv(j)] = calculateslope(data_out(:,1,j),...
            derivatives.nDeriv,derivatives.slopeorder,derivatives.slopedist,...
            plotlevel,depth_out,species{j},layercounts);
        data_out(:,2:1+derivatives.nDeriv,j)=slope;
    end
end