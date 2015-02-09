function layerdistTot = convolvelayerdist(layerdist0,layerdist1,zerolimit)

%% layerdistTot = convolvelayerdist(layerdist0,layerdist1)
% Convolve two layer number probability distributions to produce the full 
% accumulated probability distribution.
% layerdist0: layer number distribution at beginning of batch, 
% layerdist1: distribution from beginning of batch to current data point.
% All distributions have the following format:
% layerdist(:,1): Layer number
% layerdist(:,2): Corresponding probability

% Copyright (C) 2015  Mai Winstrup
% 2014-10-16 18:04: First independent script

%% Summarize probability distributions:
layerdist0 = compactdist(layerdist0,zerolimit);
layerdist1 = compactdist(layerdist1,zerolimit);

%% Maximum and minimum layer number:
yrmin = layerdist0(1,1)+layerdist1(1,1);
yrmax = layerdist0(end,1)+layerdist1(end,1); 

%% Convolve the two probability distributions:
% Initialize:
layerdistTot(:,1) = yrmin:yrmax;
layerdistTot(:,2) = 0;

for yr = yrmin:yrmax;
    % Number of layers at start of batch are layerdist0(:,1).
    % The number of annual layers in current batch must then be:
    yr_batch=yr-layerdist0(:,1);
    % Not all of these are feasible:
    mask = yr_batch>=layerdist1(1,1) & yr_batch<=layerdist1(end,1);

    % Summing the contributions:
    layerdistTot(yr-yrmin+1,2) = ...
        sum(layerdist0(mask,2).*layerdist1(yr_batch(mask)-layerdist1(1,1)+1,2));
end

%% Making sure the probabilities sum up to 1:
layerdistTot(:,2)=layerdistTot(:,2)/sum(layerdistTot(:,2));