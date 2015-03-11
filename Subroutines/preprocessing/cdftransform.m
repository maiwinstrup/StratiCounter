function ynew = cdftransform(x,y,L,plotlevel)

%% ynew = cdftransform(x,y,L,plotlevel)
% Performing a CDF-transform ("Cummulative Distribution Function") on data.
% The CDF is defined the following way: Each data point is indexed in 
% ascending order, and given a value corresponding to its index. In this 
% way, high peaks in data profile become less prominent. 
% Unless L is equal to the total data length, the transform is conducted as 
% a running algorithm over a distance of L. If L is larger than the total 
% data length, the entire data array is transformed simultaneously.
% L is given in meters.
% Copyright (C) 2015  Mai Winstrup

%% Default value of L, or if L is set as empty:
% L is then set equal to the full length of data set:
% Total data length (in meters): 
Ltot = max(x)-min(x);
if nargin==2||isempty(L); L = Ltot; end
if nargin<4; plotlevel = 0; end

%% Perform CDF-transform:
y = y(:);
ynew = nan(size(y));
if L>=Ltot;
    ynew = fcdf(y);
else
    % Same, but performed as a running mean.
    % Data around edges are found afterwards.
    istart = find(x>=x(1)+L/2,1,'first');
    iend = find(x<=x(end)-L/2,1,'last');
    for i = istart:iend
        % Perform transform for segment around data point i:
        mask = x>=x(i)-L/2 & x<=x(i)+L/2;
        segment = fcdf(y(mask));
        % Only keep the yvalue corresponding to data point i of the 
        % transformed segment:
        xmask = x(mask);
        ynew(i) = segment(xmask==x(i));
    end
    % At ends: Replaced by data of same magnitude:
    Mval = max(ynew);
    segment = y(1:istart-1);
    ynew(1:istart-1)=(segment-min(segment))/range(segment)*Mval;
    segment = y(iend+1:end);
    ynew(iend+1:end)=(segment-min(segment))/range(segment)*Mval;
end

% And finally scaling the data values to be between 0 and 1:
ynew=(ynew-min(ynew))/max(ynew);

%% Comparing the resulting data to original data profile: 
if plotlevel>1
    figure;
    subplot(2,1,1)
    plot(x,y,'-k')
    xlim([x(1) x(end)])
    title('Original data','fontweight','bold')
    subplot(2,1,2)
    plot(x,ynew)
    title('CDF transformed data','fontweight','bold')
    xlim([x(1) x(end)])
     
    % And plotting the sorted data themselves:
    figure;
    plot(sort(y))
    title('Sorted data values','fontweight','bold')
    xlabel('Index')
    ylabel('Data value')
end
end

%% Embedded cdf-transform function:
function ynew = fcdf(y)
    % Sort data according to their value:
    [~,index]=sort(y);
    index = index';
    
    % Giving data values according to index (i.e. index(1) is given a data
    % value of 1 etc):
    index(2,:)=1:length(y);

    % Now sorting according to index:
    ynew = sortrows(index',1);
    ynew = ynew(:,2)';
    
    % NaNs are kept:
    ynew(isnan(y))=nan;
end