function [logalpha_bar,logbeta,eta2,eta2_bar,gamma2,layerpos,...
    layerpos_issues,logPobs] = forwardbackward(logb,logpd,layer0_pos,...
    T,d,dmax,D,nLayerMax,tiepoints,plotlevel)

%% [logalpha_bar,logbeta,eta2,eta2_bar,gamma2,layerpos,layerpos_issues,...
%    logPobs] = forwardbackward(logb,logpd,layer0_pos,T,d,dmax,D,...
%    nLayerMax,tiepoints,plotlevel)
% Maximum a posteriori (MAP) estimation of the most likely state in a given
% pixel and corresponding uncertainties, determined using the generalized
% Forward-Backward algorithm. Layer boundary positions are given as located 
% midway between datapoints belonging to two layers. 

% Variables:
% T: length of observational data
% nLayerMax: upper-bound estimate of years in observation sequence
% layer0_pos: probability for ending of previous layer [pixels before obs]
% d: possible durations
% dmax: maximum duration
% D: number of possible durations
% logpd: corresponding log-probabilities
% logb: log-likelihood of data segments to span exactly one year
% layerpos: pixel values corresponding to mid-way between two layers, i.e. 
% do not have integer values. 

% Copyright (C) 2015  Mai Winstrup

%% Default is no plotting:
if nargin<10; plotlevel = 0; end

%% Calculate forwards and backwards variables:
% Forward message pass:
[logalpha,logalpha_bar] = forward(T,nLayerMax,layer0_pos,d,dmax,D,logpd,logb);
% Information on the previous layer position is included (if provided). 
% The previous layer ends before start of the current data batch.
% The initialization part of alpha (from before start of data batch) is 
% removed, while the last part, corresponding to t > T, is calculated and
% kept.

% Backward message pass:
logbeta = backward(T,d,dmax,logpd,logb,tiepoints,nLayerMax);
% The last part, corresponding to t > T, is kept.

% Plotting:
if plotlevel>2
    % Alpha_bar:
    figure;
    subplot(2,1,1)
    hax_alpha = gca;
    plot(hax_alpha, logalpha_bar)
    hold on
    ylabel('log \it{\alpha_t(j)}','fontweight','bold')
    xlim([1 T+dmax])
    
    % Beta: 
    subplot(2,1,2)
    hax_beta = gca;
    plot(hax_beta,logbeta)
    hold on
    ylabel('log \it{\beta_t}','fontweight','bold')
    xlim([1 T+dmax])
end

%% Posterior probability that we end a given layer at time t *and*
%% observe the given observation sequence:
% eta(t,j,d) = P(S([t-d+1:t])=j, o(1:T) | theta)
% Implemented in log-space to prevent underflow.
if isempty(tiepoints)
    logeta = logalpha+repmat(logbeta,[1 nLayerMax D]);
else
    logeta = logalpha+repmat(logbeta,[1 1 D]);
end

% Marginal probability distribution of eta over d:
% eta_bar(t,j) = sum_d(eta(t,j,d)) 
%              = P(S(t])=j, o(1:T) | theta)
logeta_bar = -ones(T+dmax,nLayerMax)*inf;
for t=1:T+dmax
    for j=1:nLayerMax
        logeta_bar(t,j)=logsumexp(logeta(t,j,:));
    end
end

%% Probability of entire observation sequence given the model parameters:
% P(o(1:T) | theta)
logPobs = logsumexp(logeta_bar(1:dmax,1));
% I.e. the joint probability of observing the observation sequence and
% ending layer 1 in the very first part of the data record (1:dmax).
% The latter has a probability equal to 1, thus this is equal to the 
% probability of the observation sequence. 
% If nLayerMax is large, we may encounter layer 1 again in the last part 
% of the data series. 

%% Posterior probability that we end a given layer at time t:
% eta2(t,j,d) = P(S([t-d+1:t])=j | obs(1,T),theta)
% eta2_bar(t,j) = P(S(t])=j | obs(1,T),theta)
eta2 = exp(logeta-logPobs);
eta2_bar = exp(logeta_bar-logPobs);

% Plotting:
if plotlevel>1
    % Probability of ending any layer at time t:
    eta2_all = sum(eta2_bar,2);
    
    figure;
    plot(eta2_bar)
    hax_eta = gca;
    hold on
    plot(eta2_all,'--k','linewidth',0.1)
    ylabel('\it{\eta_t}','fontweight','normal','fontsize',10)
    
    [X,Y]=ndgrid(1:T+dmax,1:nLayerMax);
    figure;   
    pcolor(X,Y,eta2_bar)
    shading flat
    xlabel('Pixel')
    ylabel('State')
    title('\eta_t(j)','fontweight','bold')
end

%% Probability of being in a given layer at time t:
% gamma2(t,j) = P(S(t)=j | obs(1:T),theta)
gamma2 = zeros(T,nLayerMax);

% Initialization: 
% Pixel 1 is part of layer 1:
gamma2(1,1)=1;
for t = 2:T
    gamma2(t,1)=gamma2(t-1,1)-eta2_bar(t-1,1)+eta2_bar(t-1,nLayerMax);
    % We never end layer 0 within data sequence, but we may end layer
    % nLayerMax. 
end

% Evaluated recursively for t>1 and j>1 (gamma2(1,j) is zero for j>1)
for t = 2:T
    gamma2(t,2:nLayerMax)=gamma2(t-1,2:nLayerMax)-eta2_bar(t-1,2:nLayerMax)+...
        eta2_bar(t-1,1:nLayerMax-1);
end
% Each value of gamma is calculated with a minimum precision in the order
% of 10^-16 (=eps(max(gamma))=eps(1)), and with decreasing precision in
% the tail due to the recursive evaluation. However, we are only interested
% in relatively large numbers of gamma anyway:
gamma2(gamma2<10^-5)=0;
gamma2(gamma2>1)=1;

if plotlevel>=2
    figure;
    plot(gamma2)
    hax_gamma = gca;
    hold on
    ylabel('\gamma_t','fontweight','bold')
    xlabel('\it{t}')
    
    [X,Y]=ndgrid(1:T,1:nLayerMax);
    figure;
    pcolor(X,Y,gamma2)
    xlabel('Pixel')
    ylabel('State')
    shading flat
    title('\gamma_t(j)','fontweight','bold')
end

%% The "most likely" layer boundaries (according to the MAP criterion):
[~, modevalue] = max(gamma2,[],2);
[layerpos,layerpos_issues] = findlayerposfrommode(modevalue);

% Plotting these layer boundaries onto previous plots:
if plotlevel>2
    for name = [hax_alpha, hax_beta]
        for j = 1:length(layerpos)
            plot(name,layerpos(j,1)*[1 1],[min(logbeta(:)) 0],'-k')
        end
        for j = 1:length(layerpos_issues)
            plot(name,layerpos_issues(j,1)*[1 1],[min(logbeta(:)) 0],'--r')
        end
    end
end
if plotlevel>1
    for name = [hax_eta, hax_gamma]
        hold on
        for j = 1:length(layerpos)
            plot(name,layerpos(j,1)*[1 1],[0 1],'-k')
        end
        for j = 1:length(layerpos_issues)
            plot(name,layerpos_issues(j,1)*[1 1],[0 1],'--r')
        end
    end
end
end

%% Forward variable:
function [logalpha, logalpha_bar] = ...
    forward(T,nLayerMax,layer0_pos,d,dmax,D,logpd,logb)

%% [logalpha, logalpha_bar] = forward(T,nLayerMax,layer0_pos,d,dmax,D,logpd,logb)
% Calculating the forward-variable alpha:
% log(alpha(t,j,d)) = log(P(S([t-d+1:t])=j, o(1:t) | theta))
% where theta here signifies all model parameters.
% And alpha_bar:
% log(alpha_bar(t,j)) = log(sum_d(alpha(t,j,d)))
%                = log(P(S(t])=j, o(1:t) | theta))
% I.e. alpha_bar(t,j) is the joint probability that at time t we are ending
% state j *and* that we observe the observation sequence o(1:t).
% Equations are implemented in log-space to prevent underflow.

% Variables: 
% T: length of observation sequence
% nLayerMax: upper-bound estimate of total number of years in observation sequence
% layer0_pos: probability for termination pixel of previous layer ("layer 0")
% d: possible durations
% dmax: maximum duration
% D: number of possible durations
% logpd: log-probability distribution of possible durations
% logb(t,d): log-likelihood of o(t-d+1:t) being an annual layer 

%% Initialization:
nLayerMax = nLayerMax+1; % Including data for layer 0 (index 1)
logalpha = -ones(T+2*dmax,nLayerMax,D)*Inf; % log(0)=-inf
logalpha_bar = -ones(T+2*dmax,nLayerMax)*Inf;

%% Alpha probabilities for layer 0:
logp0 = log(layer0_pos);
for i = 1:dmax % Time for ending previous layer
    logalpha(i,1,:)=logp0(i)+logpd;
    % Previous layer thickness is assumed distributed according to p(d)
    % This term is not used in subsequent summations. 
end
logalpha_bar(1:dmax,1)=logp0;

%% Recursive evaluation of alpha:
for t = 1:T+dmax
    for j = 2:nLayerMax
        % Previous alpha values:
        if j == 2
            % Layer 1: 
            % Transitions are allowed from layer 0 and layer nLayerMax 
            % (used in case of too small estimate of nLayerMax)
            if sum(isfinite(logalpha_bar(t+dmax-d,nLayerMax)))==0
                % No probabilities of nLayerMax occuring.
                logalphaPrev = logalpha_bar(t+dmax-d,j-1);
            else
                logalphaPrev = nan(D,1);
                for id = 1:D
                    logalphaPrev(id) = logsumexp([logalpha_bar(t+dmax-d(id),j-1);...
                        logalpha_bar(t+dmax-d(id),nLayerMax)]);
                end
            end
        else
            % Remaining layers:
            logalphaPrev = logalpha_bar(t+dmax-d,j-1);
        end
        
        % Set b-terms corresponding to layers starting after end of (and
        % thus being fully outside) data record to -inf. In this way, these 
        % terms amount to 0 in the sum.
        if t<=T
            logb2 = logb(t,:);
        else
            bstart = t-d+1;
            mask = bstart<=T; % Layer should start within data record
            logb2 = ones(D,1)*-inf; 
            logb2(mask)=logb(t,mask);        
        end
            
        % Calculate new alpha values:
        x = logalphaPrev+logb2(:)+logpd(:);
        logalpha(t+dmax,j,:)=x;
        logalpha_bar(t+dmax,j) = logsumexp(x);
    end
end

%% Removing initialization part of alpha and alpha_bar:
logalpha = logalpha(dmax+1:end,2:end,:);
logalpha_bar = logalpha_bar(dmax+1:end,2:end);
end

%% Backward variable:
function logbeta = backward(T,d,dmax,logpd,logb,tiepoints,nLayerMax)

%% logbeta = backward(T,d,dmax,logpd,logb,tiepoints,nLayerMax)
% Calculating the backward-variable beta:
% log(beta(t,j)) = log(sum_d(beta(t,j,d)))
%           = log(P(o(t+1:T) | S(t])=j, theta))
% where theta here signifies all model parameters.
% I.e. beta(t,j) is the probability of the given observation sequence 
% o(t+1:T), when assuming that at time t we are ending state j. Or in 
% other words - the likelihood of ending a layer at time t, given the 
% observation sequence o(t+1:T). 
% In case of no tiepoints, the boundary conditions are independent of j, 
% and thus also beta is independent on j (and d). 
% Fixpoints are introduced in such a way, that layer J (=nLayerMax) ends at 
% or after T, i.e. pixel T is part of layer J.
% Equations are implemented in log-space to prevent underflow.

% Variables: 
% T: length of observation sequence
% d: possible durations
% dmax: maximum duration
% logpd: log-probability distribution of durations
% logb(t,d): log-likelihood of o(t-d+1:t) being an annual layer 
% tiepoints: depth of possible tiepoints
% nLayerMax: maximum number of layers in batch

%% Rearrange logb-array:
% logb2(istart,iend): log-likelihood of o(istart:iend) being an annual layer 
logb2 = nan(T+2*dmax,T+2*dmax);
for tend = 1:T+dmax
    iend = tend+dmax;
    istart = iend-d+1;
    logb2(istart,iend) = logb(tend,:);
end

%% Evaluation of beta: 
% If no tiepoints:
if isempty(tiepoints)
    % Initialization:
    logbeta = zeros(T+dmax,1); % log(1)=0
    % Following (Yu, 2010, p. 219) all observations for all states are
    % taken to be equally likely for t>=T. See notes for a statistical
    % justification.
    
    % Recursive evaluation of beta for t<T:
    for t = T-1:-1:1
        x = logbeta(t+d)+logb2(t+dmax+1,t+dmax+d)'+logpd';
        logbeta(t)=logsumexp(x);
    end
    
else
    %% If tiepoints are given:
    % Initialization:
    % Layer Nmax is ending at or after T, and no other layers are allowed 
    % to end here. 
    logbeta = -ones(T+dmax,nLayerMax)*inf; % exp(-inf)=0; no probabilities
    logbeta(T:T+dmax,nLayerMax)=0; % exp(0) = 1;
    % This is according to Yu (2010, p. 219).
    
    % Recursive evaluation of beta for t<T:
    for t = T-1:-1:1
        for j = 1:nLayerMax-1
            x = logbeta(t+d,j+1)+logb2(t+dmax+1,t+dmax+d)'+logpd';
            logbeta(t,j)=logsumexp(x);
        end
    end
end
end

function y = logsumexp(x)
%% y=logsumexp(x)
% Calculate an (slightly approximate) value of y=logsumexp(x) while 
% preventing underflow in the calculations. 

%% Calculate logsumexp:
x = x(isfinite(x));
if isempty(x)
    y = -inf;
else
    % To prevent underflow:
    m = max(x);
    y = m+log(sum(exp(x-m)));
    % In case of underflow even with this manipulation, only the "large"
    % terms in the calculation are kept. The remaining terms are considered
    % insignificant anyway.(The hereby inflicted numerical error can
    % actually be decreased by normalizing to a larger constant).
end
end