function [data1,hfig] = preprocessdata(data0,depth,preprocsteps,plotlevel,...
    species,layercounts)

%% [data1,hfig] = processdata(data0,depth,preprocsteps,plotlevel,species,...
%%   layercounts)
% Preprocessing a (single) data series to obtain the best identifiable
% annual layer signal in the data. 
% The following steps may be taken (in any order): 
% Non-distance dependent transformations:
% - No processing ([])
% - Interpolation over nans ('interpNaNs')
% - Log-transformation ('log')
% - Box-Cox transformation ('boxcox',[],lambda, (alpha)])

% - Subtract constant ('minusconst',const)
% - Subtract mean ('minusmean')

% Distance-dependent transformations:
% Window length (Lwindow) is measured in meters. 
% - Normalization using quantiles ('quantile', Lwindow)
% - Normalization using min-max ('minmax', Lwindow)
% - Normalization using standard deviation ('zscore',Lwindow)
% - Using CDF-transform ('cdftransform',Lwindow)


% - Subtract baseline ('minusbaseline',[Lwindow, (quantile)])
% - Subtract smooth curve calculated using running average ('minussmooth', Lwindow)
% - Smooth using running average ('smooth',Lwindow)


% Copyright (C) 2015  Mai Winstrup

%% Set default values:
% No layercounts given, and do not display any plots. 
if nargin < 4; plotlevel=0; species=[]; layercounts=[]; end
if nargin < 5; species=[]; layercounts=[]; end
if nargin < 6; layercounts=[]; end

%% Number of processing steps:
nSteps = size(preprocsteps,1);

%% Plot original data:
hfig = gobjects(1);
if plotlevel>0
    % Depth interval shown in plot:
    [dstart_fig,dend_fig,depth_fig,data_fig] = ...
        depthrangefigure(depth,data0,layercounts);
    % Layer counts within interval:
    layercounts = layercounts(layercounts(:,1)>=dstart_fig&layercounts(:,1)<=dend_fig,:);
    % Plot original data:
    nSubfig = nSteps+1;
    hfig = plotrawdata(depth_fig,data_fig,species,nSubfig,layercounts);
    xlim([dstart_fig dend_fig])
end

%% Perform preprocessing af data in the desired order:
% Starting with the raw data:
data = data0;

for iStep = 1:nSteps
    preproctype = preprocsteps{iStep,1};
    if size(preprocsteps,2)==1
        procdist = [];
    else
        procdist = preprocsteps{iStep,2};
    end
    
    switch preproctype
        case []
            data1 = data;
              
        case 'interpNaNs'
            mask = isfinite(data);
            data1 = interp1(depth(mask),data(mask),depth);
            
        case 'log'
            % If necessary, add minimum value to prevent data from becoming 
            % negative:
            minvalue = min(data0);
            if minvalue < 0
                minvalue = floor(minvalue);
                data = data+abs(minvalue);
                disp(['Value added to ' species ' before taking log: ' num2str(abs(minvalue))])
            end
            data1 = log(data);
            
        case 'boxcox'
            lambda = procdist(1); % power parameter
            % Default for the shift parameter alpha is zero:
            if length(procdist)==1; alpha = 0;
            else alpha = procdist(2);
            end
            data1 = boxcoxtransform(data,lambda,alpha,0); % 2014-07-16 09:24
            
        case {'quantile','minmax','zscore'}
            data1 = normalizedata(depth,data,procdist,preproctype,0); % 2014-07-16 11:04
            
        case 'cdftransform'
            data1 = cdftransform(depth,data,procdist,0); % 2014-07-16 11:39

        case 'minusconst'
            data1 = data-procdist;
            
        case 'minusmean'
            data1 = data-nanmean(data0);
            
        case 'minusbaseline'
            quantiles = [0.05];
            baseline = findbaseline(depth,data0,procdist(1),quantiles,0); % 2014-07-16 11:03
            data1 = data-baseline;
         
        case 'minussmooth'
            data1 = data-smoothdata(depth,data0,procdist); % 2014-07-16 13:58

        case 'smooth'
            data1 = smoothdata(depth,data0,procdist); % 2014-07-16 13:58
    end
  
    %% Plot processed data:
    if plotlevel>0
        plotpreprocdata(hfig,nSubfig,iStep,depth,data1,preproctype,...
            procdist,layercounts,dstart_fig,dend_fig)
    end

    %% Update data:
    data0 = data1;
end
end

%% Plotting subroutines:
function hfig = plotrawdata(depth,data,species,nSubfig,layercounts)

%% plotrawdata(depth,data,species,nSubfig,layercounts):
% Plotting the original data. 

%% Plot data:
hfig(1)=figure;
subplot(nSubfig,1,1)
plot(depth,data,'-b')
hold on
plotlayercounts(layercounts,data)
title(['Before processing: ' species],'fontweight','bold','interpreter','none')
end

function plotpreprocdata(hfig,nSubfig,iStep,depth,data,proctype,procdist,...
    layercounts,dstart_fig,dend_fig)
%% plotpreprocdata(hfig,nSubfig,iStep,depth,data,proctype,procdist,...
%    layercounts,dstart_fig,dend_fig)
% This function plots a section of the processed data. 

%% Select data in depth interval:
mask = depth>=dstart_fig & depth<=dend_fig;
depth_fig = depth(mask);
data_fig = data(mask);

%% Plot data:
figure(hfig)
subplot(nSubfig,1,iStep+1)
plot(depth_fig,data_fig,'-b')
hold on
plotlayercounts(layercounts,data_fig)
title([proctype ' (' num2str(procdist) ')'],'fontweight','bold')
xlabel([dstart_fig dend_fig])
end