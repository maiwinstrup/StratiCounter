function [slope, dslope, wWhiteNoise, hfig] = calculateslope(data,order,...
    L,plotlevel,depth,species,layercounts)

%% [slope, dslope, wWhiteNoise] = calculateslope(data,order,L,plotlevel,species):
% This function calculates the slope of the intensity profile after
% Savitzky-Golay smoothing. The slope is calculated by linear regression to
% a straight line (if order=1) or parabola (if order=2), using a
% neighborhood of L pixels. The slope is then given as the differential of
% the regression curves. Different methods and neigbourhoods can be used
% for calculating the slope and the slope differential. 
% If choosing order=0, the data series is not smoothed beforehand, and
% ordinary differencing is performed.
% Derivatives are measured per pixel. 
% The derivatives are calculated per pixel, and thus their calculation 
% do not require an equidistant timescale.

% Copyright (C) 2015  Mai Winstrup
% 2014-04-21 17:07: Name updated
% 2014-05-18 19:43: Minor revisions
% 2014-07-16 14:05: showplots->plotlevel

%% Setting default values:
if size(L)==1; L=[L L]; end
if size(order)==1; order=[order 2]; end

%% Calculating the slope and slope differential for each data point:
% % Regression matrix for the slope:
% X=ones(L(1),2);
% for j =1:order(1)
%     X(:,j+1)=(-(L(1)-1)/2:(L(1)-1)/2).^j;
% end
% % And for slope differential:
% Xd=ones(L(2),2);
% for j =1:order(2)
%     Xd(:,j+1)=(-(L(2)-1)/2:(L(2)-1)/2).^j;
% end
% 
% for i = 1:N
%     % Slope:
%     % Selection of data segment:
%     px_start = max(1,i-(L(1)-1)/2);
%     px_end = min(i+(L(1)-1)/2,N);
%     datasegment = data(px_start:px_end);
%     % Selecting part of regression matrix:
%     istart = 1+px_start-i+(L(1)-1)/2;
%     Xsegment = X(istart:istart+(px_end-px_start),:);
% 
%     % Linear regression on data segment:
%     a=Xsegment\datasegment;
%     slope(i)=a(2); % Corresponding to the differential in x=0 (the mid-point)
%     
%     % Slope differential:
%     px_start = max(1,i-(L(2)-1)/2);
%     px_end = min(i+(L(2)-1)/2,N);
%     datasegment = data(px_start:px_end);
%     % Selecting part of regression matrix:
%     istart = 1+px_start-i+(L(2)-1)/2;
%     Xsegment = Xd(istart:istart+(px_end-px_start),:);
%     % Linear regression on data segment:
%     a=Xsegment\datasegment;
%     dslope(i)=2*a(3);
% end
% It is much quicker to do this by filtering as done below!

%% Savitzky-Golay smoothing:
N = length(data);
% Making sure L contains uneven numbers:
L = L-rem(L,2)+1; 
% Relative noise levels of data series:
wWhiteNoise(1)=1;

if order(1) == 0
    % Ordinary differencing:
    slope = [diff(data);nan];
    L(1)=0;
else
    % Using Savitzky-Golay smoothing before differencing:
    slope = nan(N,1);

    % Computing the 1st differential:
    [~,filtmatrix] = sgolay(order(1),L(1));
    HalfWin  = ((L(1)+1)/2)-1;
    for n = (L(1)+1)/2:length(data)-(L(1)+1)/2,
       slope(n) = dot(filtmatrix(:,2), data(n-HalfWin:n+HalfWin));
    end
    wWhiteNoise(2)=sum(filtmatrix(:,2).^2)^0.5;
end

if order(2) == 0
    dslope = [nan; diff(slope)];
    L(2)=0;
else
    % Computing the 2nd differential:
    [~,filtmatrix] = sgolay(order(2),L(2));
    HalfWin  = ((L(2)+1)/2)-1;
    dslope = nan(N,1);
    for n = (L(2)+1)/2:length(data)-(L(2)+1)/2,
        dslope(n) = 2*dot(filtmatrix(:,3)', data(n-HalfWin:n+HalfWin))';
    end
    wWhiteNoise(3)=2*sum(filtmatrix(:,3).^2)^0.5;
end

%% Plotting resulting data profile:
hfig = gobjects(1);
if plotlevel>0
    [dstart_fig,dend_fig] = selectdepthintervalforfigure(depth,data,layercounts); % based on "data" dataseries
    mask = depth > dstart_fig & depth <= dend_fig;
    
    hfig = figure;
    subplot(3,1,1)
    plot(depth(mask),data(mask))
    hold on
    plotlayercounts(layercounts,data)
    title(['Data: ' species],'fontweight','bold')
    
    subplot(3,1,2)
    plot(depth(mask),slope(mask))
    hold on
    plotlayercounts(layercounts,slope)
    title(['Slope (order=' num2str(order(1)) ', L=' num2str(L(1)) ')'],'fontweight','bold')
    
    subplot(3,1,3)
    plot(depth(mask),dslope(mask))
    hold on
    plotlayercounts(layercounts,dslope)
    title(['Curvature (order=' num2str(order(2)) ', L=' num2str(L(2)) ')'],'fontweight','bold')
    
    for i = 1:3
        subplot(3,1,i)
        xlim([dstart_fig dend_fig])
    end
end


%% Plotting subroutine:
% function hfig = plotderivatives(depth,data_final,dstart_fig,dend_fig,Model,j,counts)
% hfig = figure;
% for i = 1:3
%     subplot(3,1,i); 
%     plot(depth,data_final(:,i))
%      %
%     hold on
%     % Plot layer positions:
%     plotlayercounts(counts,data_final(:,i)) %
% end
% subplot(3,1,1)
% title(['Data: ' Model.species{j}],'fontweight','bold'); 
% subplot(3,1,2) 
% title(['Slope (order=' num2str(Model.slopeorder(1)) ', L=' num2str(Model.slopedist(1)) ')'],'fontweight','bold')
% subplot(3,1,3)
% title(['Curvature (order=' num2str(Model.slopeorder(2)) ', L=' num2str(Model.slopedist(2)) ')'],'fontweight','bold')   
% end