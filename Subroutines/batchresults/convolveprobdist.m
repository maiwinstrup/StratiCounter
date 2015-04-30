function probdistTot = convolveprobdist(probdist0,probdist1,zerolimit)

%% probdistTot = convolvelayerprob(layerdist0,layerdist1,zerolimit)
% Convolve two layer number probability distributions to produce the full 
% accumulated probability distribution.
% probdist0: layer number distribution at beginning of batch, 
% layerdist1: distribution from beginning of batch to current data point.
% All distributions have the following format:
% probdist(:,1): Layer number
% probdist(:,2): Corresponding probability

% Copyright (C) 2015  Mai Winstrup

%% Summarize probability distributions:
probdist0 = compactprobdist(probdist0,zerolimit);
probdist1 = compactprobdist(probdist1,zerolimit);

%% Maximum and minimum layer number:
yrmin = probdist0(1,1)+probdist1(1,1);
yrmax = probdist0(end,1)+probdist1(end,1); 

%% Convolve the two probability distributions:
% Initialize:
probdistTot(:,1) = yrmin:yrmax;
probdistTot(:,2) = 0;

for yr = yrmin:yrmax;
    % Number of layers at start of batch are layerdist0(:,1).
    % The number of annual layers in current batch must then be:
    yr_batch=yr-probdist0(:,1);
    % Not all of these are feasible:
    mask = yr_batch>=probdist1(1,1) & yr_batch<=probdist1(end,1);

    % Summing the contributions:
    probdistTot(yr-yrmin+1,2) = ...
        sum(probdist0(mask,2).*probdist1(yr_batch(mask)-probdist1(1,1)+1,2));
end

%% Making sure the probabilities sum up to 1:
probdistTot(:,2)=probdistTot(:,2)/sum(probdistTot(:,2));