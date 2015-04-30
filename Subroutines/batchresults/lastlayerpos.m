function postau = lastlayerpos(FBprob,tau,d,zerolimit,tstart,plotlevel)

%% postau = lastlayerpos(FBprob,tau,d,zerolimit,tstart,plotlevel)
% Calculating the probabilities corresponding to where ends the last layer 
% before the layer in pixel tau (i.e. "layer j0" of next batch).
% Copyright (C) 2015  Mai Winstrup

%% Duration parameters:
dmin = min(d);
dmax = max(d);
D = length(d);

%% Length of data batch:
batchLength = size(FBprob.gamma,1);

%% Assume layer 1 to be the most likely layer in pixel tau:
[~,j1] = max(FBprob.gamma(tau,:)); % [layer index]

%% Probabilities for where layer j0 = j1-1 ends, while knowing that layer 
% j1 covers pixel tau:
% P(S_t]=j0|S_tau=j1) = P(S_t]=j0,S_tau=j1)/P(S_tau=j1)
%                     = P(S_[t+1:tau=j1)/P(S_tau=j1)
%                     = sum_d(for t+d>=tau){P(S_[t+1:t+d]=j1)/P(S_tau=j1)}
%                     = sum_d(for t+d>=tau){eta(t+d,j1,d)}/gamma(tau,j1)

postau = zeros(dmax,1); % [pixel relative to tau]
for t = tau-dmax:tau-1
    % Summing probabilities corresponding to various values of d:
    prob=0;
    for id = max(1,(tau-t-dmin+1)):D
        prob = prob+FBprob.eta(t+d(id),j1,id);
    end
    postau(t+dmax-tau+1)=prob/FBprob.gamma(tau,j1);
end

%% Normalize probability distribution:
% First non-zero entry in postau:
px = find(postau>zerolimit,1,'first');
postau = postau(px:end); % postau is now measured in pixels before tau (tau 
% not included), assuming next batch of data to begin in pixel tau.

% Normalizing to ensure that the probability sums up to 1:
postau = postau/sum(postau);

%% Plot:
if plotlevel > 1
    plotlastlayerpos(FBprob.eta_bar,tau,postau)
    if ~isempty(tstart) % We're dealing with the actual tau for batch
        xlim([tstart batchLength])
        legend({'sum_j(\eta_{j})','\tau','pos(\tau)'})
    end
end
end

function plotlastlayerpos(eta_bar,px,pospx)
% Plot probabilities for ending position of the last layer relative to the
% one in pixel px, and compare to the total probability of ending any
% layer in this area. 

%% Probability of ending any layer at given pixel:
eta_all = sum(eta_bar,2);

%% Plot:
figure;
hline(1)=plot(eta_all,'.-b');
hold on
hline(2)=plot(px,eta_all(px),'.r','markersize',15);
hline(3)=plot((px-length(pospx)-1):px,[0;pospx;0],'-r');
xlim([px-length(pospx)-10 px+10])

xlabel('Pixel')
ylabel('Probability')
title('Probability of ending any layer j in given pixel')
legend(hline,{'sum_j(\eta_{j})',['Section boundary (px= ' num2str(px) ')'],'pos(px)'})
end