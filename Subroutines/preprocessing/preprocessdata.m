function [data1,hfig] = preprocessdata(data0,depth,preprocsteps,plotlevel,...
    species,layercounts)

%% [data1,hfig] = preprocessdata(data0,depth,preprocsteps,plotlevel,species,...
%%   layercounts)
% Preprocessing a (single) data series to obtain the best identifiable
% annual layer signal in the data. 
% The following steps may be taken (in any order): 
% Non-distance dependent transformations:
% - No processing ([])
% - Interpolation over nans ('interpNaNs')
% - Subtract mean ('minusmean')
% - Log-transformation ('log')
% - Box-Cox transformation ('boxcox',[],[lambda, (alpha)])

% Distance-dependent transformations:
% Window length (Lwindow) is measured in meters. 
% - Normalization using standard deviation ('zscore',Lwindow)
% - Normalization using min-max ('minmax', Lwindow)
% - Normalization using quantiles ('quantile', Lwindow)
% - Using CDF-transform ('cdftransform',Lwindow)
% - Subtract baseline ('minusbaseline',Lwindow,(quantile))
% - Subtract smooth curve calculated using running average ('minussmooth', Lwindow)

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
    layercounts = layercounts(layercounts(:,1)>=dstart_fig&...
        layercounts(:,1)<=dend_fig,:);
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
    procdist = [];
    procval = [];
    if size(preprocsteps,2)>=2; procdist = preprocsteps{iStep,2}; end
    if size(preprocsteps,2)==3; procval = preprocsteps{iStep,3}; end
    
    switch preproctype
        case []
            data1 = data;
              
        case 'interpNaNs'
            mask = isfinite(data);
            data1 = interp1(depth(mask),data(mask),depth);
            
        case 'minusmean'
            data1 = data-nanmean(data);
            
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
            % Power parameter (lambda):
            lambda = procval(1); 
            % Default for the shift parameter (alpha) is zero:
            if length(procval)==1; alpha = 0;
            else alpha = procval(2);
            end
            data1 = boxcoxtransform(data,lambda,alpha,0); 
            
        case {'zscore','minmax','quantile'}
            data1 = normalizedata(depth,data,procdist,preproctype,0); 
            
        case 'cdftransform'
            data1 = cdftransform(depth,data,procdist,0); 

        case 'minusbaseline'
            quantile = 0.05; % baseline corresponds to the 5% quantile
            baseline = findbaseline(depth,data0,procdist(1),quantile,0);
            data1 = data-baseline;
         
        case 'minussmooth'
            data1 = data-smoothdata(depth,data,procdist); 
    end
  
    %% Plot processed data:
    if plotlevel>0
        plotpreprocdata(hfig,nSubfig,iStep,depth,data1,preproctype,...
            procdist,procval,layercounts,dstart_fig,dend_fig)
    end

    %% Update data:
    data0 = data1;
end
end

%% Plotting subroutines:
function hfig = plotrawdata(depth,data,species,nSubfig,layercounts)

%% plotrawdata(depth,data,species,nSubfig,layercounts):
% Plot the original data. 

%% Plot data:
hfig(1)=figure;
subplot(nSubfig,1,1)
plot(depth,data,'-b')
hold on
plotlayercounts(layercounts,data)
title(['Before processing: ' species],'fontweight','bold','interpreter','none')
end

function plotpreprocdata(hfig,nSubfig,iStep,depth,data,preproctype,procdist,...
    procval,layercounts,dstart_fig,dend_fig)
%% plotpreprocdata(hfig,nSubfig,iStep,depth,data,proctype,procdist,...
%    procval,layercounts,dstart_fig,dend_fig)
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
% Title:
text = preproctype;
if ~isempty(procdist); text = [text '(' num2str(procdist) ')']; end
if ~isempty(procval); text = [text '(' num2str(procval) ')']; end
title(text,'fontweight','bold','interpreter','none')
xlim([dstart_fig dend_fig])
end