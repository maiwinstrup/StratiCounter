function [slope, derivnoise, hfig] = calculateslope(data,nDeriv,...
    slopeorder,slopedist,plotlevel,depth,species,layercounts)

%% [slope, derivnoise, hfig] = calculateslope(data,nDeriv,slopeorder,...
%    slopedist,plotlevel,depth,species,layercounts)
% This function calculates derivatives of the intensity profile using 
% Savitzky-Golay smoothing. The slope is calculated by linear regression to
% a straight line (if slopeorder=1) or parabola (if slopeorder=2), using a
% neighborhood of slopedist pixels. Slope is calculated based on the 
% differential of the regression curves. As derivatives are calculated per 
% pixel, their calculation does not require data on equidistant depth scale.
% Copyright (C) 2015  Mai Winstrup

%% Set default values:
% Default is to calculate 1st derivative using linear polynomium, 2nd 
% derivatives using 2nd degree polynomium etc.
if isempty(slopeorder); slopeorder=1:nDeriv; end

% If only one value given for slope distance, this value should be used 
% for all calculations of derivatives:
if size(slopedist)==1; slopedist=ones(1,nDeriv)*slopedist; end
% Make sure that slopedist contains uneven numbers:
slopedist = slopedist-rem(slopedist,2)+1; 

%% Calculating the slope and slope differential for each data point:
% % Regression matrix for the slope:
% X=ones(slopedist(1),2);
% for j =1:slopeorder(1)
%     X(:,j+1)=(-(slopedist(1)-1)/2:(slopedist(1)-1)/2).^j;
% end
% % And for slope differential:
% Xd=ones(slopedist(2),2);
% for j =1:slopeorder(2)
%     Xd(:,j+1)=(-(slopedist(2)-1)/2:(slopedist(2)-1)/2).^j;
% end
% 
% N = length(data)
% for i = 1:N
%     % Slope:
%     % Selection of data segment:
%     px_start = max(1,i-(slopedist(1)-1)/2);
%     px_end = min(i+(slopedist(1)-1)/2,N);
%     datasegment = data(px_start:px_end);
%     % Selecting part of regression matrix:
%     istart = 1+px_start-i+(slopedist(1)-1)/2;
%     Xsegment = X(istart:istart+(px_end-px_start),:);
% 
%     % Linear regression on data segment:
%     a=Xsegment\datasegment;
%     Slope corresponding to the differential in x=0 (the mid-point)
%     slope(i)=a(2); 
%     
%     % Slope differential:
%     px_start = max(1,i-(slopedist(2)-1)/2);
%     px_end = min(i+(slopedist(2)-1)/2,N);
%     datasegment = data(px_start:px_end);
%     % Selecting part of regression matrix:
%     istart = 1+px_start-i+(slopedist(2)-1)/2;
%     Xsegment = Xd(istart:istart+(px_end-px_start),:);
%     % Linear regression on data segment:
%     a=Xsegment\datasegment;
%     dslope(i)=2*a(3);
% end
% It is much quicker to calculate this by filtering as done below!

%% Savitzky-Golay smoothing:
% Relative noise levels of data series is defined as equal to one:
derivnoise(1)=1;

% Calculate derivatives: 
N = length(data);
slope = nan(N,nDeriv);

% Using Savitzky-Golay smoothing before differencing. 
% Compute the 1st differential:
if nDeriv >= 1
    [~,filtmatrix] = sgolay(slopeorder(1),slopedist(1));
    halfwin  = ((slopedist(1)+1)/2)-1;
    for i = (slopedist(1)+1)/2:N-(slopedist(1)+1)/2,
        slope(i,1) = dot(filtmatrix(:,2), data(i-halfwin:i+halfwin));
    end
    % Relative noiselevel of data series:
    derivnoise(2)=sum(filtmatrix(:,2).^2)^0.5;
end

% Calculate 2nd order differential:
if nDeriv >= 2
    [~,filtmatrix] = sgolay(slopeorder(2),slopedist(2));
    halfwin  = ((slopedist(2)+1)/2)-1;
    for i = (slopedist(2)+1)/2:N-(slopedist(2)+1)/2,
        slope(i,2) = 2*dot(filtmatrix(:,3)', data(i-halfwin:i+halfwin))';
    end
    % Relative noise level:
    derivnoise(3)=2*sum(filtmatrix(:,3).^2)^0.5;
end

%% Plot data and derivatives:
hfig = gobjects(1);
if plotlevel>0&&nDeriv>0
    % Select depth interval for figure:    
    [dstart_fig,dend_fig] = depthrangefigure(depth,data,layercounts); 
    
    % Layers and data in this interval:
    layercounts = layercounts(layercounts(:,1)>=dstart_fig&layercounts(:,1)<=dend_fig,:);
    mask = depth>=dstart_fig & depth<=dend_fig;
    
    % Plot original data:
    hfig = figure;
    subplot(nDeriv+1,1,1)
    plot(depth(mask),data(mask))
    hold on
    plotlayercounts(layercounts,data(mask))
    title(['Data: ' species],'fontweight','bold','interpreter','none')
    xlim([dstart_fig dend_fig])
    
    % Name of derivative data series: 
    name = {'Slope','Curvature'};
    for k = 1:nDeriv
        subplot(nDeriv+1,1,k+1)
        plot(depth(mask),slope(mask,k))
        hold on
        plotlayercounts(layercounts,slope(mask,k))
        title([name{k} ' (order=' num2str(slopeorder(k)) ...
            ', dist=' num2str(slopedist(k)) ')'],'fontweight','bold')
        xlim([dstart_fig dend_fig])
    end
end