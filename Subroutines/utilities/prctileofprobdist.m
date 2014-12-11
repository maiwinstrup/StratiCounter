function q = prctileofprobdist(x,p,prc)

%% q = prctileofprobdist(x,p,prc)
% Calculating percentiles (prc, in values<1) of a probability distribution. 
% x must be a sorted column of x-values with associated probabilities p. 
% If p is a matrix, each column of the matrix are considered probability
% distributions for the x-values.

% Mai Winstrup
% 2014-08-08 10:52

%% Normalize probabilities: 
p = p./repmat(sum(p,1),size(p,1),1);

%% CDF-matrix of probability distributions (extended):
pj = [zeros(1,size(p,2)); cumsum(p,1); ones(1,size(p,2))];

% Corresponding x-values:
xj = [x(1)-1; x; x(end)+1];
     
%% Calculating percentiles:
q = nan(length(prc),size(p,2));
for i = 1:size(p,2)
    for j = 1:length(prc)
        q(j,i) = xj(find(pj(:,i)>prc(j),1,'first'));
    end
end