function hfig = showtimescale(timescale,timescale1yr,manualcounts,dbatch,...
    Model,filename)

%% hfig = showtimescale(timescale,timescale1yr,manualcounts,dbatch,...
%   Model,filename)
% Plot the resulting timescale, and its associated uncertainties, as 
% function of depth. If manual counts exist, these are plotted along with 
% their estimated uncertainties. In this way, the automated layer-counted 
% timescale can easily be compared to the manual counts. 

% Copyright (C) 2015  Mai Winstrup

%% Initialize figure: 
hfig=figure;
clf
set(hfig, 'paperunits', 'centimeters','papertype', 'A4', ...
   'paperorientation', 'portrait','paperposition',[1 1 11 8])
set(hfig, 'DefaultAxesFontSize', 9);
set(hfig, 'DefaultAxesUnits','normalized')

%% Plot manual layer counts (if exist): 
if ~isempty(manualcounts)
    % Remove manual counts from outside plotting range: 
    mask = manualcounts(:,1)>=Model.dstart & ...
        manualcounts(:,1)<=timescale(end,1);
    manualcounts = manualcounts(mask,:);
end

% Counts in current interval?
if ~isempty(manualcounts)
    % Convert layer counts into stepped data files:
    [layerpos_manual, age_manual] = ...
        steplayerpos(manualcounts(:,1),manualcounts(:,2),Model);
    
    % Stepped version of accumulated uncertainty from beginning of section:
    unc_manual = [manualcounts(1,4) manualcounts(1:end,4)'; ...
        manualcounts(1,4) manualcounts(1:end,4)'];
    unc_manual = unc_manual(:);
    unc_manual = unc_manual - unc_manual(1);
else
    layerpos_manual = nan;
    age_manual = nan;
    unc_manual = nan;
end

% Plot uncertainty band:
color_manual  = [1 1 1]*0.5;
alpha_manual = 0.7;
xvalues = [layerpos_manual+Model.dx/10; layerpos_manual(end:-1:1)];
% dx is added only for plotting purposes
yvalues_unc = [age_manual-unc_manual-0.01; age_manual(end:-1:1)+unc_manual(end:-1:1)+0.01];
% 0.01 is added/subtracted only for plotting purposes
% Number of points used to create line (forth and back):
N = min(length(xvalues),2000); % Maximum ~2000 points
dx = round(length(xvalues)/N);
fill(xvalues(1:dx:end),yvalues_unc(1:dx:end),color_manual,'edgecolor',...
    color_manual,'facealpha',alpha_manual,'edgealpha', alpha_manual)
hold on
% Plot the most likely ages:
hline(1)=plot(layerpos_manual(1:dx:end),age_manual(1:dx:end),'-k','linewidth',1);
legendname{1} = Model.nameManualCounts(1:end-4); % remove fileextension

%% Timescale results from Forward-Backward algorithm: 
% Mode and confidence intervals of age-distribution at each datapoint.
% Confidence interval bands: 
color  = [1 0.3 0.3];
alpha = 0.5;
% Number of confidence intervals:
nConf = (size(timescale,2)-2)/2;

% Confidence bands are plotted on top of each other:
for i = 1:nConf
    unc_lowerbound = timescale(:,2+i);
    unc_upperbound = timescale(:,2+2*nConf-i+1);    
    xvalues = [timescale(:,1)+Model.dx/10; timescale(end:-1:1,1)];
    yvalues = [unc_lowerbound-0.01; unc_upperbound(end:-1:1)+0.01];

    % Reduce array to only include those points where yvalues change:
    mask = diff(yvalues)~=0;
    mask = [true;mask]|[mask;true]; % Selecting points on both sides of boundary
    xvalues = xvalues(mask);
    yvalues = yvalues(mask);

    % Number of points used to create line (forth and back):
    N = min(length(xvalues),2000); 
    dx = round(length(xvalues)/N);
    fill(xvalues(1:dx:end),yvalues(1:dx:end),color,'edgecolor',color,...
        'facealpha',alpha,'edgealpha', alpha)    
end

% Plot most likely solution, i.e. mode of the distribution:
hline(2)=plot(timescale(1:dx:end,1),timescale(1:dx:end,2),...
    'color',[0.9 0 0],'linewidth',1.5);
legendname{2} = 'Forward-Backward';

%% Results from timescale1yr:
% Convert layer positions into stepped files:
[layerpos_1yr, age_1yr] = steplayerpos(timescale1yr(:,1),timescale1yr(:,2),Model);
% And plot these:
hline(3)=plot(layerpos_1yr(1:dx:end),age_1yr(1:dx:end),'-',...
    'color',color*0.6,'linewidth',1.2);
legendname{3} = [Model.icecore '_timescale'];

%% Mark the depths of batch boundaries:
% If tiepoints: Batch boundaries are the location of tiepoints.
ylimit = get(gca,'ylim');
K = length(legendname);
if ~isempty(dbatch)
    hline(K+1) = plot(dbatch(1)*[1 1],ylimit,'--','color',[1 1 1]*0.8);
    legendname{K+1} = 'Batch divisions';
    for i = 2:length(dbatch)
        plot(dbatch(i)*[1 1],ylimit,'--','color',[1 1 1]*0.8)    
    end
end

%% Mark the depths of marker horizons, if any:
K = length(legendname);
if ~isempty(Model.dMarker)
    legendname{K+1} = 'Marker horizons';
    for iMarkerSet = 1:length(Model.dMarker)
        for i = 1:length(Model.dMarker{iMarkerSet})
            dMarker = Model.dMarker{iMarkerSet}(i);
            
            % Size of marker on plot:
            if strcmp(Model.ageUnitOut,'AD')
                indexmin = size(timescale,2); 
                indexmax = 3;
            else
                indexmin = 3;
                indexmax = size(timescale,2);
            end
            
            if isempty(manualcounts)
                ymin = interp1(timescale(:,1),timescale(:,indexmin),dMarker);
                ymax = interp1(timescale(:,1),timescale(:,indexmax),dMarker);
            else
                ymin = nanmin(interp1(manualcounts(:,1),...
                    manualcounts(:,2)-manualcounts(:,4),dMarker,'linear',nan),...
                    interp1(timescale(:,1),timescale(:,indexmin),dMarker));
                ymax = nanmax(interp1(manualcounts(:,1),...
                    manualcounts(:,2)+manualcounts(:,4),dMarker,'linear',nan),...
                    interp1(timescale(:,1),timescale(:,indexmax),dMarker));
            end
            
            % Plot line:
            hline(K+1) = plot(dMarker*[1 1],[ymin-5 ymax+5],'-','color',[1 1 1]*0.4);
        end
    end
end

%% Add legend and axis labels:
% Labels and limits:
h_ax = gca; 
xlabel(h_ax,'Depth [m]','fontweight','bold')
ylabel(h_ax,['Age (' Model.ageUnitOut ')'],'fontweight','bold')
title([Model.icecore ' layer counts'],'fontweight','bold','interpreter','none')
xlim([Model.dstart timescale(end,1)])
box on
% Legend:
legend(hline,legendname,'location','Best','fontsize',7,'interpreter','none')
legend('boxoff')

%% Save figure:
if ~isempty(filename)
    print(filename,'-djpeg','-r400')    
end
end

%% Embedded stepping function:
function [dlayer_step, age_step] = steplayerpos(dlayer,age,Model)
% For easy plotting: Step the arrays of layer depth positions and 
% corresponding ages, according to selected age unit "ageUnitOut", 

% Layer depths:
dlayer = dlayer(:);
dlayer_step = [dlayer(1)-1 dlayer(1:end)'; dlayer(1:end)' dlayer(end)+1];
dlayer_step = dlayer_step(:);

% Ages:
age = age(:);
switch Model.ageUnitOut
    case 'AD'
        age_step = [age(1:end)' age(end)-1; age(1:end)' age(end)-1];
    otherwise
        age_step = [age(1)-1 age(1:end)'; age(1)-1 age(1:end)'];
end
age_step = age_step(:);
end