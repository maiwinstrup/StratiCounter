function  [Layerpos,FBprob,ExpVal,logPobs,d,pd,logb_tot,nLayerMax_new] = ...
    layerdetection(data,Template,Layerpar,Layer0,nLayerMax,Model,plotlevel)

%% [Layerpos, FBprob, ExpVal, logPobs, d, pd, logb_tot, bweight] = 
%   layerdetection(data,Template,Layerpar,Layer0,nLayerMax,Model,plotlevel)
% This script calculates the most likely layer for each pixel of a data 
% batch, and the corresponding uncertainties. These probabilities are 
% computed using the Forward-Backward algorithm. Layer template and layer 
% parameters employed are given (“Template“, “Layerpar”), as well as initial 
% conditions (“Layer0”) and possible constraints in form of tiepoints 
% (Model.tiepoints).
% An additional output of the algorithm is the estimation of several 
% quantities (ExpVal), later to be used for calculating the likelihood of 
% the applied layer model and parameters when conditioned on the observed 
% data.

% Copyright (C) 2015  Mai Winstrup
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.

%% Layer thickness probabilities: 
% Layer durations (thicknesses) are assumed to follow either a lognormal 
% distribution or gamma distribution.
% Definition: Pixel 1 covers the interval ]0.5,1.5]. 
% The continous probability distributions are thus discretized according to 
% P(i)=int(p(x)dx,i-0.5,i+0.5), where P: discrete probs, p: continuous probs.
% In this way, a given pixel will belong to the layer that covers the most
% (more than half) of the pixel.

% Possible state durations (d) and corresponding probabilities:
if strcmp(Model.durationDist,'logn')
    % Upper bound for d:
    meanLambda = exp(Layerpar.my+Layerpar.sigma^2/2); %[m]
    dmax = 10*meanLambda/Model.dx; %[pixel] 
    % Cumulative probability distribution:
    pcum = logncdf((0.5:1:dmax+0.5)*Model.dx,Layerpar.my,Layerpar.sigma);

elseif strcmp(Model.durationDist,'gamma')
    % Upper bound for d:
    disp('gamma distribution not supported yet');
    meanLambda = Layerpar.shape*Layerpar.scale; %[m]
    dmax = 10*meanLambda/Model.dx; %[pixel] 
    % Cumulative probability distribution:
    pcum = gamcdf((0.5:1:dmax+0.5)*Model.dx,Layerpar.shape,Layerpar.scale);
end 

% Removing tails of the distribution, equally much in both ends.
% Approximately Model.tailPercent are removed at each end. 
dmin = max(1,find(pcum>Model.tailPrc*10^-2,1,'first')-1); 
dmax = find(pcum<1-Model.tailPrc*10^-2,1,'last'); 
% In this way, equally many should be over and under-estimated. Observe
% that any given small value of lambda removed from consideration is far
% more likely than a larger one. 

% Error message if pcum does not cover an appropriate range of the
% probability distribution:
if dmax == length(pcum)
    disp(['OBS: Spread of duration distribution is large. ' ...
        'Significant part of the tail is removed.'])
end

% Normalizing duration probabilities:
% To ensure the probabilities sum up to 100%, the remaining tails are
% either distributed evenly or spliced onto the ends.
if strcmp(Model.tailType,'evenly') % Recommended
    % Even distribution of remaining probabilies onto all durations.
    % Duration probabilities:
    pd = pcum(dmin+1:dmax+1) - pcum(dmin:dmax);
    pd = pd/sum(pd);    
elseif strcmp(Model.tailType,'spliced')
    % Additional probabilities are spliced onto the ends of the
    % distribution:
    pd = [pcum(dmin+1:dmax) 1] - [0 pcum(dmin+1:dmax)];
end 

% Possible durations:
d = dmin:1:dmax;
% Total number of options:
D = length(d);

% Plot duration distribution:
if plotlevel>2
    figure;
    plot(d*Model.dx*100,pd)
    title('Duration probability function','fontweight','bold')
    xlabel('d [cm]')
    ylabel('p(d)')
end

%% Initial conditions: 
% Probability distribution for ending of previous layer:
% Vector of length dmax providing the probabilities of ending the previous 
% layer in these pixels prior to start of the current batch. 
layer0_pos = zeros(1,dmax);
istart = length(Layer0.pos)-dmax+1;
layer0_pos(max(1,-istart+2):end)=Layer0.pos(max(1,istart):end); 
% The value of dmax may have changed: In that case, the probability of 
% layer boundaries in the very first pixels - removed above - are added to 
% the probabilities corresponding to the first pixel.
layer0_pos(1)=sum(Layer0.pos(1:istart));
% Normalize probabilities:
layer0_pos = layer0_pos/sum(layer0_pos);

%% Exhaustive calculation of the layer likelihood, b:
% b(t1,t2) = P(obs(max(t1,1):min(t2,T)) | S([t1:t2])=j)
% is calculated for all values of t1 and t1+min(d)-1<=t2=<t1+max(d)-1
% All layers are assumed to be an outcome of the same process, and thus "j"
% stands for any layer.
% Also calculated are various expectation values later required for
% calculating an updated estimate of the layer parameters.
% NaNs in data sequence are allowed.

% Introducing index notation: 
% Data is appended with dmax NaNs in the beginning and/or end of data
% series. Indices (i) and data number (t) are related as: i = t+dmax

% Length of data record:
T = size(data,1);

% There is no need to calculate b-values for values of t1 corresponding to 
% indices for which the previous layer (layer0) is known to not yet have 
% ended (i.e. layer0_pos equals zero). The first possible starting point 
% for the first layer (layer1) in the current batch is:
i0 = find(layer0_pos>0,1,'first')+1; 

if strcmp(Model.bcalc,'BLR')
    % Using Bayesian Linear Regression to a generic layer template:
    [logb, ExpVal] = likelihoodBLR(data,Template,Layerpar,d,i0,Model); 
    % logb(:,:,j) contains probabilities for the various chemical species (j).  
    % The calculated probability values are log-probability densities and 
    % may therefore have values above 0 (i.e. b>1).
    % ExpVal provides a range of expectation values later to be used for 
    % estimating a new set of layer parameters.

    % Treatment of areas with no-data in data records: 
    % Replacing by the mean probability value, i.e. we assumed average
    % probability (for given species) of layer boundaries in these sections:
    nodatamask = logb==0; % Corresponding to a probability of 1
    datamask = ~nodatamask & isfinite(logb);
    for j = 1:Model.nSpecies
        speciesmask = false(size(logb));
        speciesmask(:,:,j)=true;
        % Mean probability value:
        meanvalue = mean(logb(speciesmask & datamask)); 
        if isfinite(meanvalue) % If no data exists in interval: keep value of 1
            logb(nodatamask & speciesmask) = meanvalue;
        end
    end
    
    % The all-species layer likelihood is the weighted sum of the layer 
    % likelihood contributions from the various species:
    wSpecies(1,1,:) = Model.wSpecies(:);
    wSpeciesMatrix = repmat(wSpecies,size(logb(:,:,1)));
    logb_tot = sum(logb.*wSpeciesMatrix,3);

    % Accounting for the multitude of data series: 
    % The general assumption is that the relative probabilities of shape vs. 
    % layer thickness should be weighted as: 
    % bweight = 1/nEff,
    % with the effective number of data series (nEff) being dependent on 
    % their weighting values, and the number of derivative data series: 
    nEff = (Model.derivatives.nDeriv+1)*sum(Model.wSpecies);
    
    % However, we may want to use another weighting (or adjust it). 
    % If so, the weighting can be changed by choosing a value of 
    % Model.bweight different from 1. 
    bweight = Model.bweight/nEff;
    % A value of Model.bweight>1: The dependency on layer shapes are 
    % enhanced relative to layer durations. 
    
    % Weighted all-species layer likelihood:
    logb_tot = logb_tot*bweight;
else
    disp('This calculation of log(b) is not implemented')
end

% Plot layer likelihood and expectation values:
if plotlevel>=2
    % log(b):
    hfig = figure;
    hax = nan(1,Model.nSpecies+2);
    for j = 1:Model.nSpecies
        hax(j) = subplot(1,Model.nSpecies+1,j);
        convertandplot(logb(:,:,j),d,dmax,T)
        title(Model.species{j},'fontweight','bold')
    end
    hax(j+1)=subplot(1,Model.nSpecies+1,Model.nSpecies+1);
    convertandplot(logb_tot,d,dmax,T)
    title('log(b_{tot})*w','fontweight','bold')
    suptitle('log(b)')
        
    % Compute and plot log(b*p(d)):
    logbpd = nan(T+dmax,D);
    for id = 1:D
        logbpd(:,id) = logb_tot(:,id)+log(pd(id));
    end
    figure;
    hax(j+2)=gca;
    convertandplot(logbpd,d,dmax,T)
    title('log(b_{tot}^{w}*p(d))','fontweight','bold')
    
    % Plot estimated layer parameters:
    if plotlevel > 2
        plotlayerpar(ExpVal,Layerpar,Model,d,dmax,D,T);
    end
end

%% The Forward-Backward algorithm:
% Using the Forward-Backward algorithm to obtain a Maximum-A-Posteriori (MAP) 
% estimate of the most likely state in a given pixel and corresponding 
% uncertainties, along with the probability of the data sequence given the 
% current model parameters. 
[FBprob.logalpha_bar,FBprob.logbeta,FBprob.eta,FBprob.eta_bar, ...
    FBprob.gamma,Layerpos.fb,Layerpos.fb_issues,logPobs] = ...
    forwardbackward(logb_tot,log(pd),layer0_pos,T,d,dmax,D,nLayerMax,...
    Model.tiepoints,plotlevel);

% Plot resulting layer boundaries onto figures:
if plotlevel>=2
    addlayerstofigs(hax,Layerpos,Model,dmax)
end

%% Is the possible number of layers in batch larger than anticipated? 
% If so, increase value of nLayerMax in next iteration.
if max(FBprob.gamma(:,end))>10^-3 % If the probability of being in 
    % layer j=nMax at any point t is larger than this value.
    % Number of times that all layers numbers have been circled through:
    nCircles = round(sum(FBprob.eta_bar(:,end)));
    % Largest possible end layer value:
    [~,nLayerEnd] = find(FBprob.gamma(end,:)>0,1,'last');    
    nLayerMax_new = ceil(1.1*(nLayerMax*nCircles+nLayerEnd));
    disp(['nLayerMax increased:' num2str(nLayerMax_new)])
else
    nLayerMax_new = nLayerMax;
end
end

%% Embedded plotting subroutines:
function convertandplot(matrix,d,dmax,T)
%% convertandplot(matrix,d,dmax,T)
% Convert from matrix(tend,id) to matrix(istart,iend), and plot results. 

%% Convert matrix(tend,id) to matrix2(istart,iend):
matrix2 = nan(T+2*dmax,T+2*dmax);
for tend = 1:T+dmax
    iend = tend+dmax;
    istart = iend-d+1;
    matrix2(istart,iend) = matrix(tend,:);
end

%% Plot matrix2:
hold on
pcolor(matrix2')
shading flat
grid on
xlim([1 T+2*dmax])
ylim([1 T+2*dmax])

% Indicate area corresponding to data series:
plot([dmax dmax],[1 T+2*dmax],'-k')
plot([T+dmax T+dmax],[1 T+2*dmax],'-k')
plot([1 T+2*dmax],[dmax dmax],'-k')
plot([1 T+2*dmax],[T+dmax T+dmax],'-k')

% Line corresponding to a layer of length 1:
plot(1:T+2*dmax,1:T+2*dmax,'-k')
xlabel('t_{start}')
ylabel('t_{end}')
colorbar
set(gca, 'xtick', dmax:100:T+dmax-1,'xticklabel', 1:100:T,...
    'ytick', dmax:100:T+dmax-1,'yticklabel', 1:100:T)
end

function hfig = plotlayerpar(ExpVal,Layerpar,Model,d,dmax,D,T)
%% hfig = plotlayerpar(ExpVal,Layerpar,Model,d,dmax,D,T) 
% Plot estimated most likely layer parameters.      
for j = 1:Model.nSpecies
    Er = reshape([ExpVal(:,:,j).r]',T+dmax,D,Model.order);
    hfig(j) = figure;
    for ipar = 1:Model.order
        subplot(1,Model.order+1,ipar)
        convertandplot(Er(:,:,ipar)+Layerpar.par(ipar),d,dmax,T)
        title(['\phi_' num2str(ipar) '^0 + E(r_' num2str(ipar) ')'],'fontweight','bold')
    end
        
    % Residuals: 
    E_nvar =nan(T+dmax,D);
    for id = 1:D
        for tend = 1:T+dmax
            if ~isempty(ExpVal(tend,id,j).X)
                E_nvar(tend,id)=sum((ExpVal(tend,id,j).oXr-ExpVal(tend,id,j).X*Layerpar.par(:,j)).^2);
            end
        end
    end
    subplot(1,Model.order+1,Model.order+1)
    convertandplot(E_nvar,d,dmax,T)
    title('e^2,  e=o_j-X(\phi^0 + E(r))','fontweight','bold')
    suptitle(Model.species{j})
end
end

function addlayerstofigs(hax,Layerpos,Model,dmax)
%% addlayerstofigs(hax,Layerpos,Model,dmax)
% Add layer boundaries to figures.

%% Plotting the resulting layer boundaries onto the calculated log(b):
for k = 1:Model.nSpecies+2
    axes(hax(k))
    hold on
    for j = 1:length(Layerpos.fb(:,1))-1
        % Forward-Backward layer boundaries:
        plot(Layerpos.fb(j,1)+dmax,Layerpos.fb(j+1,1)+dmax, '.k','markersize',7)
    end
end
end