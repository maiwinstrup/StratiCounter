function [data1,hfig,dstart_fig,dend_fig] = processdata(data,species,preprocess,Model,depth,counts,plotlevel)

%% [data1,hfig,dstart_fig,dend_fig] = processdata(data,species,preprocess,Model,depth,counts,plotlevel)
% Preprocessing data to obtain the best identifiable annual layer signal
% in the data series. The following steps may be taken (in any order): 
% - No preprocessing ('none', [])
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
% "Depth" and "counts" are only used for plotting.

% Mai Winstrup
% 2014-04-21 17:12: Filename updated
% 2014-05-13 01:22: Set default values
% 2014-05-18 20:08: Depth is nolonger part of data array
% 2014-07-15 16:56: Plotlevel introduced, new format of preprocess
% 2014-07-16 16:18: Unequidistant data allowed, downsamling and derivatives 
%                   moved to processdata.m
% 2014-08-15 17:08: output for hfig when plotlevel=0: [] -> nan
% 2014-08-17 22:00: minusmean: min->nanmean
% 2014-08-23 13:55: check of depthscale moved from script

%% Set default values:
% No counts given, do not display any plots. 
if nargin < 6; counts = []; end
if nargin < 7; plotlevel = 0; end

%% Plot original data:
if plotlevel > 0
    [hfig,nSubfig,dstart_fig,dend_fig] = plotdata_raw(depth,data,species,preprocess,Model,counts);
else
    dstart_fig = nan; 
    dend_fig = nan;
    hfig = [];
end

%% Perform preprocessing af data in the desired order:
% Starting with the raw data:
data0 = data;
% Number of processing steps:
nSteps = size(preprocess,1);

for iStep = 1:nSteps
    preproc_type = preprocess{iStep,1};
    const = preprocess{iStep,2};
    
    switch preproc_type
        case 'none'
            data1 = data0;
            
        case 'interpNaNs'
            mask = isfinite(data0);
            data1 = interp1(depth(mask),data0(mask),depth);
            
        case 'log'
            % If necessary: Add minimum value to prevent data from becoming negative
            minvalue = min(data0);
            if minvalue < 0
                minvalue = floor(minvalue);
                data0 = data0+abs(minvalue);
                disp(['Value added to ' species ' before taking log: ' num2str(abs(minvalue))])
            end
            data1 = log(data0);
            
        case 'boxcox'
            lambda = const(1); % power parameter
            % Default for the shift parameter alpha is zero:
            if length(const)==1; alpha = 0;
            else alpha = const(2);
            end
            data1 = boxcoxtransform(data0,lambda,alpha,0); % 2014-07-16 09:24
            
        case {'quantile','minmax','zscore'}
            data1 = normalizedata(depth,data0,const,preproc_type,0); % 2014-07-16 11:04
            
        case 'cdftransform'
            data1 = cdftransform(depth,data0,const,0); % 2014-07-16 11:39

        case 'minusconst'
            data1 = data0-const;
            
        case 'minusmean'
            data1 = data0-nanmean(data0);
            
        case 'minusbaseline'
            baseline = findbaseline(depth,data0,const(1),const(2),0); % 2014-07-16 11:03
            data1 = data0-baseline;
         
        case 'minussmooth'
            data1 = data0-smoothdata(depth,data0,const); % 2014-07-16 13:58

        case 'smooth'
            data1 = smoothdata(depth,data0,const); % 2014-07-16 13:58
    end
  
    %% Plot processed data:
    if plotlevel > 0
        plotdata_proc(hfig,nSubfig,iStep,depth,data1,preproc_type,const,dstart_fig,dend_fig,counts)
    end

    %% Update data:
    data0 = data1;
end
end

%% Plotting subroutines:
function [hfig,nSubfig,dstart_fig,dend_fig] = plotdata_raw(depth,data,species,...
    preprocess,Model,counts)
%% plotdata_raw(data,depth,species,preprocess,Model,counts):
% This function plots a section of the original data. The section is chosen
% to cover approximately 25 layers. 
% Mai Winstrup, 2014

%% Depth section to be shown in plot: 
dstart_fig = Model.dstart;
if ~isempty(counts)
    L = mean(diff(counts(:,1)))*25; % Approximately 25 layers   
else
    L = 1;
end
dend_fig = min(dstart_fig+L,Model.dend);
mask = depth>=dstart_fig & depth<=dend_fig;
% Depthscale for figure:
depth_fig = depth(mask);
data_fig = data(mask);

% Data should cover at least half of section:
while sum(isnan(data_fig))>0.5*length(depth_fig)
    dstart_fig = dend_fig;
    dend_fig = min(dstart_fig+L,Model.dend);
    mask = depth>=dstart_fig & depth<=dend_fig;
    depth_fig = depth(mask);
    data_fig = data(mask);
    
    % If this never is possible: Use data starting from the beginning. 
    if dend_fig==Model.dend 
        dstart_fig = Model.dstart;
        dend_fig = Model.dstart+L;
        break; 
    end
end

%% Plot data:
% Number of subfigures:
nSubfig = size(preprocess,1)+1;
hfig(1)=figure;
subplot(nSubfig,1,1)
plot(depth_fig,data_fig,'-b')
hold on
plotlayercounts(counts,data_fig)
xlim([dstart_fig dend_fig])
title(['Before preprocessing: ' species],'fontweight','bold','interpreter','none')
end

function plotdata_proc(hfig,nSubfig,iStep,depth,data,preproc_type,const,dstart_fig,dend_fig,counts)
% This function plots a section of the processed data. The depth section is
% the same as for the plot of the raw data. 
% Mai Winstrup, 2014

%% Select data in depth interval:
mask = depth>=dstart_fig & depth<=dend_fig;
depth_fig = depth(mask);
data_fig = data(mask);

%% Plot data:
figure(hfig)
subplot(nSubfig,1,iStep+1)
if strcmp(preproc_type,'interpNaNs')
    plot(depth_fig,data_fig,'--b')
else
    plot(depth_fig,data_fig,'-b')
end
hold on
plotlayercounts(counts,data_fig)
xlim([dstart_fig dend_fig])
title([preproc_type ' (' num2str(const) ')'],'fontweight','bold')
end