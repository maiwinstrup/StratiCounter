function baseline = findbaseline(x,y,L,q,plotlevel)

%% FINDBASELINE(x,y,L,q,plotlevel):
% Calculates a baseline for the data array. The baseline is defined as the
% q'th quantile of the distribution of the data values in an interval of
% length L around each data point. Multiple baselines are allowed.

% Mai Winstrup, 2011
% 2014-05-18 19:08: Minor adjustments
% 2014-07-16 11:03: showplots->plotlevel,unequidistant data allowed

%% Set default value of q to 0.5:
if isempty(q); q = 0.5; end

%% Calculating baseline:
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
    plot(y)
    hold on
    plot(baseline,'-k')
    title('Data and calculated baseline(s)','fontweight','bold')
end