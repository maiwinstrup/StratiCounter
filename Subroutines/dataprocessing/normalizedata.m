function ynew = normalizedata(x,y,L,method,plotlevel)

%% ynew = normalizedata(x,y,L,method,plotlevel)
% Running normalization of the signal contained in y(x). The normalization 
% is made over a distance of L (measured in same units as x), and can be 
% done in 3 different ways: 'minmax','zscore' and 'quantile'. NaNs in data
% series are allowed.

% Copyright (C) 2015  Mai Winstrup
% 2014-05-18 20:08: Minor adjustments
% 2014-07-16 11:04: showplots->plotlevel, allowing for unequidistant data

%% Normalization: 
N = length(y);
ynew = nan(N,1);

%% Quantile normalization: 
if strcmp(method,'quantile')
    % Normalization to constant probability distribution via quantiles
    q = 0:0.1:1;
    baseline = findbaseline(x,y,L,q,plotlevel); % 2014-07-16 11:03
    for i = 1:N
        if isfinite(y(i))
            ynew(i) = interp1q(baseline(i,:)',q',y(i));
        end
    end
    return
end

%% Min-max/zscore normalization:
for i = 1:N
    % Segment around data point i:
    mask = x>=x(i)-L/2 & x<=x(i)+L/2;
    ysegment = y(mask);
    
    switch method
        case 'minmax'
            % Min-Max normalization:
            % Normalization to between 0 and 1.
            ynew(i) = (y(i)-nanmin(ysegment))/(nanmax(ysegment)-nanmin(ysegment));
            
        case 'zscore'
            % Z-score normalization:
            % Normalization to constant mean and standard deviation
            ynew(i)=(y(i)-nanmean(ysegment))/nanstd(ysegment);
    end
end