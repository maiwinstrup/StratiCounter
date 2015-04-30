function probdist = compactprobdist(probdist,zerolimit)

%% probdist = compactprobdist(probdist,zerolimit)
% This function compacts the probability distribution with entries:
% dist(:,1): Values of distribution, dist(:,2): probabilities. 
% This is done by: 1) removing tails of distribution, incl. entries of 
% zero, and 2) renormalizing the distribution.

% Copyright (C) 2015  Mai Winstrup

%% Normalize distribution: 
probdist(:,2)=probdist(:,2)/sum(probdist(:,2));

%% Remove tails of distribution: 
% Discard a tiny amount of tails of the distribution: 
index1 = find(cumsum(probdist(:,2))>zerolimit,1,'first');
index2 = find(cumsum(probdist(:,2))>=1-zerolimit,1,'first'); 
if isempty(index2); index2 = max(probdist(:,1)); end

% The probability contained in the tails is spliced onto the ends:
probdist(index1,2)=probdist(index1,2)+sum(probdist(1:index1-1,2));
probdist(index2,2)=probdist(index2,2)+sum(probdist(index2+1:end,2));
% Tails are removed:
probdist = probdist(index1:index2,:);

%% Re-normalizing the probability distribution:
probdist(:,2)=probdist(:,2)/sum(probdist(:,2));