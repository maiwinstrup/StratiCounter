function dist = compactdist(dist,zerolimit)

%% dist = compactdist(dist,zerolimit)
% This function compacts the probability distribution with entries:
% dist(:,1): Values of distribution, dist(:,2): probabilities. 
% This is done by: 1) removing tails of distribution, incl. entries of 
% zero, and 2) renormalizing the distribution.

% Mai Winstrup,
% 2014-10-14 15:32: Added normalization before removal of tails
% 2014-10-16 13:51: zerolimit introduced as variable

%% Normalize distribution: 
dist(:,2)=dist(:,2)/sum(dist(:,2));

%% Remove tails of distribution: 
% Discard a tiny amount of tails of the distribution: 
index1 = find(cumsum(dist(:,2))>zerolimit,1,'first');
index2 = find(cumsum(dist(:,2))>=1-zerolimit,1,'first'); 
if isempty(index2); index2 = max(dist(:,1)); end

% The probability contained in the tails is spliced onto the ends:
dist(index1,2)=dist(index1,2)+sum(dist(1:index1-1,2));
dist(index2,2)=dist(index2,2)+sum(dist(index2+1:end,2));
% Tails are removed:
dist = dist(index1:index2,:);

%% Re-normalizing the probability distribution:
dist(:,2)=dist(:,2)/sum(dist(:,2));