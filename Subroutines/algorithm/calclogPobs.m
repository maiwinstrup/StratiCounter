function logPobs_new = calclogPobs(logPobs,Layerpar,Prior,Model)
% Copyright (C) 2015  Mai Winstrup

%% Value of logPobs when also accounting for priors:
if strcmp(Model.update{1},'QB')
    % Probability of found value of my: 
    logpmy = -0.5*log(2*pi*Prior.v)-(Layerpar.my-Prior.m)^2/(2*Prior.v);
else
    logpmy = 0;
end
    
if strcmp(Model.update{4},'QB')
    % Probability of found parameter values:
    logppar = nan(1,Model.nSpecies);
    for j = 1:Model.nSpecies
        L = chol(inv(Prior.invU(:,:,j)),'lower');
        v = L\(Layerpar.par(:,j)-Prior.u(:,j));
        logppar(j) = -Model.order/2*log(2*pi)-sum(log(diag(L)))-0.5*sum(v.^2);
    end
else
    logppar = 0;
end
        
% Resulting posterior probability of observing data sequence, when
% accounting for priors:
logPobs_new = logPobs + logpmy + sum(logppar);