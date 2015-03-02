function [data1,hfig] = preprocessdata(data,depth,procsteps,plotlevel,...
    species,layercounts)

%% [data1,hfig,dstart_fig,dend_fig] = processdata(data,depth,procsteps,...
%    plotlevel,species,layercounts)
% Preprocessing data (a single data series) to obtain the best identifiable
% annual layer signal in the data series. 
% The following steps may be taken (in any order): 
% - No processing ([])
% - Interpolation over nans ('interpNaNs', [])
% - Log-transformation ('log', [])
% - Box-Cox transformation ('boxcox', [lambda, (alpha)])
% - Normalization using quantiles ('quantile', Lwindow)
% - Normalization using min-max ('minmax', Lwindow)
% - Normalization using standard deviation ('zscore',Lwindow)
% - Using CDF-transform ('cdftransform',Lwindow)
% - Subtract constant ('minusconst',const)
% - Subtract mean ('minusmean')
% - Subtract baseline ('minusbaseline',[Lwindow, (quantile)])
% - Subtract smooth curve calculated using running average ('minussmooth', Lwindow)
% - Smooth using running average ('smooth',Lwindow)
% Window length (Lwindow) is measured in meters. 
% "Depth" and "layercounts" are only used for plotting.

% Copyright (C) 2015  Mai Winstrup
% 2014-04-21 17:12: Filename updated
% 2014-05-13 01:22: Set default values
% 2014-05-18 20:08: Depth is nolonger part of data array
% 2014-07-15 16:56: Plotlevel introduced, new format of proc
% 2014-07-16 16:18: Unequidistant data allowed, downsamling and derivatives 
%                   moved to processdata.m
% 2014-08-15 17:08: output for hfig when plotlevel=0: [] -> nan
% 2014-08-17 22:00: minusmean: min->nanmean
% 2014-08-23 13:55: check of depthscale moved from script

%% Set default values:
% No counts given, do not display any plots. 
if nargin < 3; plotlevel=0; species=[]; layercounts=[]; depth=[]; end

%% Plot original data:
hfig = gobjects(1);
if plotlevel > 0
    % Depth interval shown in plot:
    [dstart_fig,dend_fig,depth_fig,data_fig] = selectdepthintervalforfigure(depth,data,layercounts);

    nSubfig = size(procsteps,1)+1;
    hfig = plotrawdata(depth_fig,data_fig,species,nSubfig,layercounts,dstart_fig,dend_fig);
else
    dstart_fig = nan; 
    dend_fig = nan;
end

%% Perform preprocessing af data in the desired order:
% Starting with the raw data:
data0 = data;
% Number of processing steps:
nSteps = size(procsteps,1);

for iStep = 1:nSteps
    proctype = procsteps{iStep,1};
    if size(procsteps,2)==1
        procdist = [];
    else
        procdist = procsteps{iStep,2};
    end
    
    switch proctype
        case []
            data1 = data0;
              
        case 'interpNaNs'
            mask = isfinite(data0);
            data1 = interp1(depth(mask),data0(mask),depth);
            
        case 'log'
                % If necessary: Add minimum value to prevent data from becoming 
            % negative
            minvalue = min(data0);
            if minvalue < 0
                minvalue = floor(minvalue);
                data0 = data0+abs(minvalue);
                disp(['Value added to ' species ' before taking log: ' num2str(abs(minvalue))])
            end
            data1 = log(data0);
            
        case 'boxcox'
            lambda = procdist(1); % power parameter - og altså ikke en afstand!!
            % Default for the shift parameter alpha is zero:
            if length(procdist)==1; alpha = 0;
            else alpha = procdist(2);
            end
            data1 = boxcoxtransform(data0,lambda,alpha,0); % 2014-07-16 09:24
            
        case {'quantile','minmax','zscore'}
            data1 = normalizedata(depth,data0,procdist,proctype,0); % 2014-07-16 11:04
            
        case 'cdftransform'
            data1 = cdftransform(depth,data0,procdist,0); % 2014-07-16 11:39

        case 'minusconst'
            data1 = data0-procdist;
            
        case 'minusmean'
            data1 = data0-nanmean(data0);
            
        case 'minusbaseline'
            baseline = findbaseline(depth,data0,procdist(1),procdist(2),0); % 2014-07-16 11:03
            data1 = data0-baseline;
         
        case 'minussmooth'
            data1 = data0-smoothdata(depth,data0,procdist); % 2014-07-16 13:58

        case 'smooth'
            data1 = smoothdata(depth,data0,procdist); % 2014-07-16 13:58
    end
  
    %% Plot processed data:
    if plotlevel > 0
        plotprocdata(hfig,nSubfig,iStep,depth,data1,proctype,procdist,dstart_fig,dend_fig,layercounts)
    end

    %% Update data:
    data0 = data1;
end
end

%% Plotting subroutines:
function [hfig,dstart_fig,dend_fig] = plotrawdata(depth,data,...
    species,nSubfig,layercounts,dstart_fig,dend_fig)
%% plotrawdata(data,depth,species,procsteps,layercounts):
% This function plots a section of the original data. The section is chosen
% to cover approximately 25 layers. 
% Mai Winstrup

%% Select data in depth interval:
mask = depth>=dstart_fig & depth<=dend_fig;
depth_fig = depth(mask);
data_fig = data(mask);

%% Plot data:
% Number of subfigures:
hfig(1)=figure;
subplot(nSubfig,1,1)
plot(depth_fig,data_fig,'-b')
hold on
plotlayercounts(layercounts,data_fig)
xlim([dstart_fig dend_fig])
title(['Before processing: ' species],'fontweight','bold','interpreter','none')
end

function plotprocdata(hfig,nSubfig,iStep,depth,data,proctype,procdist,...
    dstart_fig,dend_fig,layercounts)
% This function plots a section of the processed data. The depth interval 
% is the same as for the plot of the raw data. 
% Mai Winstrup

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
xlim([dstart_fig dend_fig])
title([proctype ' (' num2str(procdist) ')'],'fontweight','bold')
end