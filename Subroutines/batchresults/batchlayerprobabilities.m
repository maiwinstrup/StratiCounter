function [ntau, ntauTotal, LayerProbDist] = ...
    batchlayerprobabilities(depth,FBprob,tau,Layer0,prctile,plotlevel)

%% [ntau, ntauTotal, LayerProbDist] = batchlayerprobabilities(depth,FBprob,...
%   tau,Layer0,prctile,plotlevel)
% Calculating the probability distribution of number of layers in current 
% batch (ntau) and total number of layers in data series up to end of 
% current batch (ntauTotal), along with various properties of the changing 
% probability distribution along batch of total number of layers in the 
% data series (LayerProbDist).

% Copyright (C) 2015  Mai Winstrup

%% Resulting number of annual layers in current data batch:
% Probability distribution of annual layer number (j) in pixel tau (i.e. in 
% the beginning of a new layer):
% P(n(tau)) = P(S(tau)=j|obs(1:T),theta) = gamma(tau,j)

% Maximum number of layers in batch (incl. start layer):
nLayerMax = size(FBprob.gamma,2);
% Layer numbers, and corresponding probabilities:
ntau = nan(nLayerMax,2);
ntau(:,1)=0:nLayerMax-1;
ntau(:,2)=FBprob.gamma(tau,:);
% Note: First layer is here counted as layer 0, such that we are in fact 
% counting "layer boundaries". This is appropriate since the batch division
% has been selected such that we start in the same layer as the previous 
% batch ends (right after a well-defined layer boundary), and this layer 
% should not be counted twice. 

% Compact the distribution (remove tails, incl. entries of zero, and 
% re-normalize the resulting distribution):
zerolimit = 10^-5;
ntau = compactprobdist(ntau,zerolimit); 

% Plot the layer distribution in pixel tau: 
if plotlevel>=2
    plotlayerprobfortau(ntau)
    title(['Layer probability for \tau=' num2str(tau) 'px'],'fontweight','bold')
end

%% Total annual layer probability distribution through batch:
% All layers (and their uncertainties) through the first batch to the 
% current one is accounted for. This is done via convolving with the 
% probabilities given in Layer0.number, providing the full layer number 
% probability distribution in the first pixel of the batch. 
% The combined timescale starts with layer number 0, i.e. we are counting 
% "layer boundaries".

% Total minimum and maximum number of layers by end of current batch:
yrmin = min(Layer0.no(:,1));
yrmax = max(Layer0.no(:,1))+ntau(end,1); 
zerolimit = 10^-6;

% Full probability distribution for each pixel: 
layerprob = zeros(tau,yrmax-yrmin+1);
for t = 1:tau
    % Resulting probability distribution:
    probdist = convolveprobdist(Layer0.no,...
        [(0:ntau(end,1))',FBprob.gamma(t,1:ntau(end,1)+1)'],zerolimit); 
    % Inset in layer probability array:
    layerprob(t,probdist(:,1)-yrmin+1) = probdist(:,2);
    % Note: Indexing in layerprob (and gamma) is layer+1, since index 0 
    % does not exist, whereas indexing in ntau etc. is according to layer
    % number in column 1, which may contain layer number 0. 
end

%% Full, accumulated, layer distribution in pixel tau:
% Layer numbers and corresponding probabilities:
ntauTotal(:,1)=yrmin:yrmax;
ntauTotal(:,2)=layerprob(tau,:);
% Compact the distribution (remove tails, renormalize etc):
zerolimit = 10^-5;
ntauTotal = compactprobdist(ntauTotal,zerolimit); 

% Plot the layer distribution in pixel tau: 
if plotlevel>=2
    plotlayerprobfortau(ntauTotal)
    title('Layer probability for \tau (all batches)','fontweight','bold')
end

%% Mode, mean, median and quantiles of annual layer number for each pixel:
% These values are based on the accumulated probability distribution that 
% accounts for layering and uncertainties in all previous batches.
% Calculated for t<=tau.

% Corresponding depthscale:
LayerProbDist.d = depth(1:tau);

% Mode of distribution:
[~, LayerProbDist.mode]=max(layerprob,[],2);
LayerProbDist.mode = LayerProbDist.mode(:)+yrmin-1; % The '-1' is due to indexing

% Calculating running mean, median and percentiles:
LayerProbDist.mean=nan(tau,1);
LayerProbDist.median=nan(tau,1);
LayerProbDist.prctile=nan(tau,length(prctile));

for t = 1:tau
    % Mean of distribution at t:
    LayerProbDist.mean(t) = sum(layerprob(t,:).*(yrmin:yrmax));
    % Median:
    LayerProbDist.median(t) = prctileofprobdist((yrmin:yrmax)',(layerprob(t,:))',0.5);
    % Percentiles: 
    LayerProbDist.prctile(t,:)=prctileofprobdist((yrmin:yrmax)',(layerprob(t,:))',prctile);
end

if plotlevel>2
    plotlayerstats(LayerProbDist)
end
end

%% Subfunctions for plotting
function plotlayerprobfortau(ntau)
%% plotlayerprobfortau(ntau)
% Plot layer number distribution (ntau) in pixel tau.

% Stepped probabilities:
[nstep, pstep] = stepit([ntau(1,1)-1;ntau(:,1);ntau(end,1)+1],[0;ntau(:,2);0]);
% Plot:
figure;
plot(nstep,pstep, '-k', 'linewidth',1)
xlabel('Layer number')
ylabel('Probability')
end

function plotlayerstats(LayerProbDist)
%% plotlayerstats(LayerProbDist)
% Plot layerdistribution along current batch, using various definitions of
% the most likely layer (mean, median, mode), and the percentiles of the
% distribution. 

figure;
% Mode:
hline(1)=plot(LayerProbDist.d,LayerProbDist.mode, '-k','linewidth',2);
hold on
% Percentiles of distribution:
for i = 1:size(LayerProbDist.prctile,2)
    hline(2)=plot(LayerProbDist.d,LayerProbDist.prctile(:,i), 'color',[1 1 1]*0.5,'linewidth',1);
end
% Median: 
hline(3)=plot(LayerProbDist.d,LayerProbDist.median(:,1), '-b');
% Mean:
hline(4)=plot(LayerProbDist.d,LayerProbDist.mean(:,1), '-g');
title('Mode, mean and median of layer number distribution','fontweight','bold')
hlegend = legend(hline,{'Mode','Percentiles','Median','Mean'});
set(hlegend,'location','best');
xlabel('Depth')
ylabel('Layer number')
end

function [xstep, ystep]=stepit(x,y)
%% [xstep, ystep]=stepit(x,y)
% A function to step data for plotting purposes
% Sune Olander Rasmussen

dx = (diff(x))';
x = [x(1)-mean(dx)/2,x(2:end)'-dx/2; x(1:end-1)'+dx/2, x(end)+mean(dx)/2];
xstep = x(:);

y = [y(:)'; y(:)'];
ystep = y(:);
end