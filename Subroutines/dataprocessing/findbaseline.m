function baseline = findbaseline(x,y,L,q,plotlevel)

%% findbaseline(x,y,L,q,plotlevel):
% Calculates a baseline for the data array. The baseline is defined as the
% q'th quantile of the distribution of the data values in an interval of
% length L around each data point. Multiple baselines are allowed.
% Copyright (C) 2015  Mai Winstrup

%% Set default values:
% No plotting:
if nargin < 5; plotlevel = 0; end
% Set default value of baseline equal to the 5% quantile;
if nargin < 4 || isempty(q); q = 0.05; end

%% Calculate baseline:
N = length(y);
baseline = nan(N,length(q));
for i=1:N
    mask = x>=x(i)-L/2 & x<=x(i)+L/2;
    ysegment = y(mask);
    baseline(i,:)=quantile(ysegment,q);
end

%% Plotting:
if plotlevel > 0
    figure;
    plot(x,y)
    hold on
    plot(x,baseline,'-k')
    title('Data and calculated baseline(s)','fontweight','bold')
end