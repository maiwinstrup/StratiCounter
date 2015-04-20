function [data_out,depth_out,derivnoise,hfigpreproc,hfigderiv] = ...
    makedatafile(data_in,depth_in,preprocsteps,derivatives,dx,dx_center,...
    plotlevel,species,layercounts)

%% [data_out,depth_out,derivnoise,hfigpreproc,hfigderiv] = 
%      makedatafile(data_in,depth_in,preprocsteps,derivatives,dx,dx_center,
%      plotlevel,species,layercounts)
% Take in a data series, which may contain multiple species, and preprocess 
% these according to the specifications in preprocsteps. Subsequently, the 
% data series may be downsampled to equidistant depthscale (unless dx=[]),
% and derivatives are calculated. "Derivnoise" is a vector providing the 
% theoretical amount of white noise in the derivative data series relative 
% to the original data profile. Preprocessed data and derivatives are 
% plotted if plotlevel>0. Figure handles are provided as output. 
% Copyright (C) 2015  Mai Winstrup

%% Set default: 
if nargin<9; layercounts=[]; end % No plotting of layer counts
if nargin<8; species=[]; end; % No titles on plot
if nargin<7; plotlevel=0; end % No plotting
if nargin<6; dx_center=0; end  % Pixels are not centered in depth
if nargin<5; dx=[]; end  % No adjustment of depth scale

% If not given, the species name array is empty:
nSpecies = size(data_in,3);
if isempty(species)
    for j = 1:nSpecies; species{j}=''; end
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

%% Interpolate depth scale:
if isempty(dx)
    depth_out = depth_in;
    data_out(:,1,:)=data1(:,1,:);
    % Initialize rest of array:
    data_out(:,2:1+derivatives.nDeriv,:)=nan;

else
    % Interpolate to equidistant depth scale: 
    % Species #1:
    [depth_out,data_out(:,1,1)] = ...
        downsampling(depth_in,data1(:,1,1),dx,dx_center);    
    % Initialize data array for remaining species (and derivatives):
    L = length(depth_out);
    data_out = cat(3,[data_out(:,1,1),nan(L,derivatives.nDeriv)], ...
        nan(L,1+derivatives.nDeriv,nSpecies-1));
    
    % Remaining species:
    for j = 2:nSpecies
        [depth_out,data_out(:,1,j)] = ...
            downsampling(depth_in,data1(:,1,j),dx,dx_center);
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

    elseif isempty(dx) && isempty(preprocsteps{j})
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