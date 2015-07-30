function ynew = normalizedata(x,y,L,method,plotlevel)

%% ynew = normalizedata(x,y,L,method,plotlevel)
% Running normalization of the signal contained in y(x). The normalization 
% is made over a distance of L (measured in same units as x), and can be 
% done in 3 different ways: 'zscore', 'minmax' and 'quantile'. NaNs in data
% series are allowed.
% Copyright (C) 2015  Mai Winstrup

%% Set default: No plotting
if nargin < 5; plotlevel = 0; end

%% Initialization: 
ynew = nan(size(y)); % Same shape as input file

%% Min-max/zscore normalization:
for i = 1:length(y)
    % Segment around data point i:
    mask = x>=x(i)-L/2 & x<=x(i)+L/2;
    ysegment = y(mask);
    
    switch method
        case 'zscore'
            % Z-score normalization:
            % Normalization to constant mean and standard deviation
            ynew(i)=(y(i)-nanmean(ysegment))/nanstd(ysegment);

        case 'minmax'
            % Min-Max normalization:
            % Normalization to between 0 and 1.
            ynew(i) = (y(i)-nanmin(ysegment))/(nanmax(ysegment)-nanmin(ysegment));            
       
    end
end

%% Quantile normalization: 
if strcmp(method,'quantile')
    % Normalization to constant probability distribution via quantiles
    q = 0:0.2:1;
    baseline = findbaseline(x,y,L,q,plotlevel); 
    for i = 1:length(y)
        if isfinite(y(i))
            ynew(i) = interp1q(baseline(i,:)',q',y(i));
        end
    end
end