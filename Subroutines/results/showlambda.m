function hfig = showlambda(lambda,Layerpos,timescale1yr,manualcounts,Model,filename)

%% hfig = showlambda(lambda,Layerpos,timescale1yr,manualcounts,Model,filename)
% Plot mean annual layer thicknesses as averaged over a given equidistant 
% depth intervals. Layerthicknesses derived from both manual and 
% automatically layer-counted timescales are plotted. The resulting plot is 
% saved under "filename". 

% Copyright (C) 2015  Mai Winstrup
% 2014-08-13 22:30: Using new version of lambda array
% 2014-08-14 11:13: Adjusted multiplyfactor
% 2014-08-15 11:12: Removed .jpeg from filename

%% If lambda is empty (or very small)
% hfig is empty, and return to main. 
if size(lambda,1)<=2; hfig = nan; return; end

%% Initialize figure: 
hfig = figure;
set(hfig, 'paperunits', 'centimeters','papertype', 'A4', ...
   'paperorientation', 'portrait','paperposition',[1 1 15 10])
set(hfig, 'DefaultAxesFontSize', 10, 'DefaultAxesUnits','normalized')
h = 0.8; w = 0.8; 
dy0 = (1-h)/2;
dx0 = (1-w)/2;
hax = axes('position',[dx0 dy0 w h],'nextplot','add','yaxislocation','left','xaxislocation','bottom');

%% Select unit for layer thickness: 
meanLambda = mean(lambda(:,3));
if meanLambda < 10^-3;
    lambdaunit = 'm';
    multiplyfactor = 1;
elseif meanLambda <10^-2
    lambdaunit = 'mm';
    multiplyfactor = 10^3;
elseif meanLambda < 1
    lambdaunit = 'cm';
    multiplyfactor = 10^2;
else
    lambdaunit = 'm';
    multiplyfactor = 1;
end

%% Mean layer thickness of manually counted layers:
% If no manual counts exist: Add NaNs
if isempty(manualcounts)
    manualcounts(1,1:4)=nan;
end

% Section boundaries (these may be different than for auto-counts in 
% beginning and end of data)
lambda_man = lambda(:,1:2);
lambda_man(lambda_man<manualcounts(1,1)-meanLambda) = manualcounts(1,1);
lambda_man(lambda_man>manualcounts(end,1)+meanLambda) = manualcounts(end,1);

% Number of years, a.k.a. layer boundaries, within each section:
nManual = nan(size(lambda,1),1); 
% Maximum Counting Error:
mce = nan(size(lambda,1),1);

for i = 1:size(lambda,1)
    mask = manualcounts(:,1)>=lambda_man(i,1) & manualcounts(:,1)<lambda_man(i,2);
    if sum(mask)==0
        nManual(i)=nan;
        mce(i)=nan;
    else
        % Number of years:
        % Uncertain years are counted as half-years
        nManual(i) = sum(mask)-sum(manualcounts(mask,3))*0.5; 
        % Maximum counting error over interval:
        mce(i) = sum(manualcounts(mask,3))*0.5;
    end
end   

% Layer thicknesses:
Lman = lambda_man(:,2)-lambda_man(:,1);
lambdaManual(:,1) = Lman./nManual;
lambdaManual(:,2) = Lman./(nManual-mce); % Minimum value
lambdaManual(:,3) = Lman./(nManual+mce); % Maximum value
% If no layer boundaries, the mean layer thickness is NaN

% Plot onto figure:
dplot = [lambda_man(:,1)'; lambda_man(:,2)'];
dplot = dplot(:);
% Uncertainty bands:
xvalues = [dplot; dplot(end:-1:1)];
lambdaManualMin = [lambdaManual(:,2)'; lambdaManual(:,2)'];
lambdaManualMax = [lambdaManual(:,3)'; lambdaManual(:,3)'];
lambdaManualMax = lambdaManualMax(:);
yvalues = [lambdaManualMin(:); lambdaManualMax(end:-1:1)];
color_gicc  = [1 1 1]*0.5;
alpha_gicc = 0.7;
fill(xvalues(:),yvalues(:)*multiplyfactor,color_gicc,'edgecolor',color_gicc,'facealpha',alpha_gicc,'edgealpha', alpha_gicc)
hold on
% Most likely layer thickness:
lambdaManualML = [lambdaManual(:,1)'; lambdaManual(:,1)']; 
% Plot in appropriate unit:
hline(1)=plot(dplot,lambdaManualML(:)*multiplyfactor,'-k','linewidth',2);
legendtext{1} = 'Manual';

%% Mean layer thicknesses from Forward-Backward algorithm as run for each 
% section individually:
lambdaFBmode = [lambda(:,3)'; lambda(:,3)']; 
dplot = [lambda(:,1)'; lambda(:,2)'];
dplot = dplot(:);
hline(2)=plot(dplot,lambdaFBmode(:)*multiplyfactor,'-r','linewidth',2);
hold on

% Uncertainty bands for each confidence interval: 
nConf = length(Model.prctile)/2;
for i = 1:nConf
    lambdaFBmin = [lambda(:,3+i)'; lambda(:,3+i)'];
    lambdaFBmax = [lambda(:,4+2*nConf-i)'; lambda(:,4+2*nConf-i)'];
    lambdaFBmax = lambdaFBmax(:);
    xvalues = [dplot; dplot(end:-1:1)];
    yvalues = [lambdaFBmin(:); lambdaFBmax(end:-1:1)];
    color = [1 0.3 0.3];
    alpha = 0.5;
    fill(xvalues(:),yvalues(:)*multiplyfactor,color,'edgecolor',color,'facealpha',alpha,'edgealpha',0)
end
legendtext{2} = 'Forward-Backward';

%% Mean layer thicknesses from output when calculated based on layer 
% boundaries:
% Using the improved Forward-Backward boundaries from timescale1yr:
nAuto = nan(size(lambda,1),1);
for i = 1:size(lambda,1)
    mask=timescale1yr(:,1)>=lambda(i,1)&timescale1yr(:,1)<lambda(i,2);
    nAuto(i)=sum(mask);
end
L = lambda(:,2)-lambda(:,1);
lambdaAutoBoundaries = L./nAuto;
lambdaAutoBoundaries = [lambdaAutoBoundaries'; lambdaAutoBoundaries'];   
hline(3) = plot(dplot,lambdaAutoBoundaries(:)*multiplyfactor,'--','color',color*0.6,'linewidth',1.2);
legendtext{3} = 'Layers in timescale1yr';

%% Layer thicknesses from Viterbi output:
% if ~isempty(Layerpos.viterbi)
%     nViterbi = nan(size(lambda,1),1);
%     for i = 1:size(lambda,1)
%         mask=Layerpos.viterbi>=lambda(i,1)&Layerpos.viterbi<lambda(i,2);
%         nViterbi(i) = sum(mask);
%     end
%     lambdaViterbi = L./nViterbi;
%     lambdaViterbi = [lambdaViterbi'; lambdaViterbi'];   
%     hline(4)=plot(dplot,lambdaViterbi(:)*multiplyfactor,'-m','linewidth',2);
%     legendtext{4} = 'Viterbi';
% end

%% Tiepoints as grey bars:
ymin = min(lambda(:,3))*multiplyfactor;
ymax = max(lambda(:,end))*multiplyfactor;
K = length(legendtext);
if ~isempty(Model.tiepoints)
    % Find tiepoints within selected depth interval:
    mask = Model.tiepoints(:,1)>=lambda(1,1)&Model.tiepoints(:,1)<=lambda(end,2);
    dtiepoints = Model.tiepoints(mask,1);
    if ~isempty(dtiepoints)
        legendtext{K+1} = 'Tiepoints';
        for i = 1:length(dtiepoints)
            hline(K+1)=plot(dtiepoints(i)*[1 1],[ymin ymax],'--','color',[1 1 1]*0.8,'linewidth',4);
        end
    end
end

%% Mark the depths of marker horizons, if any:
K = length(legendtext);
if ~isempty(Model.dMarker)
    for iMarkerSet = 1:length(Model.dMarker)
        % Find horizons within selected depth interval:
        mask = Model.dMarker{iMarkerSet}>=lambda(1,1)&Model.dMarker{iMarkerSet}<=lambda(end,2);
        dMarker = Model.dMarker{iMarkerSet}(mask);
        if ~isempty(dMarker)
            legendtext{K+1} = 'Marker horizons';
            for i = 1:length(dMarker)
                hline(K+1) = plot(dMarker(i)*[1 1],[ymin ymax],'-','color',[1 1 1]*0.4*i/length(dMarker));
            end
        end
    end
end

%% Add legend, title and labels:
% Legend:
hleg = legend(hline,legendtext);
set(hleg,'box','off')
% Label:
ylabel(['Mean layer thickness [' lambdaunit ']']) 
xlabel('Depth [m]','fontweight','bold')
% Title:
dxText = mode(L);
if dxText < 1
    title(['Calculated over ' num2str(dxText*100) 'cm sections'],'fontweight','bold')
else
    title(['Calculated over ' num2str(dxText) 'm sections'],'fontweight','bold')
end

%% Add timescale on top (based on auto counts)
% Only to be done if covering more than 500 years.
tinterval = range(timescale1yr(:,2));
if tinterval>500
    % Spacing of ticks: 
    dt = roundsignificant(tinterval/5,1);
    ttick_start = ceil(min(timescale1yr(:,2))/dt)*dt;
    ttick_end = floor(max(timescale1yr(:,2))/dt)*dt;
    ttick = ttick_start:dt:ttick_end;
    xticklabel = cell(length(ttick),1);
    for i = 1:length(ttick) 
        xticklabel{i} = [num2str(ttick(i)) ' ' Model.ageUnitOut];
    end
    % Depth of the tick marks:
    dtick = interp1(timescale1yr(:,2),timescale1yr(:,1),ttick);
    % Plot in dummy-axes:
    hax0 = axes('position',[dx0 dy0 w h],'nextplot','add',...
        'xaxislocation','top','ytick',[],'yticklabel',[],'color','none');
    set(hax0,'xlim',get(hax,'xlim'))
    % Ensure increasing numbers: 
    [dtick,index]=sort(dtick,'ascend');
    xticklabel = xticklabel(index);
    set(hax0,'xtick',dtick,'xticklabel',xticklabel)
end

%% Save figure:
if ~isempty(filename)
    print(filename,'-djpeg','-r400')    
end