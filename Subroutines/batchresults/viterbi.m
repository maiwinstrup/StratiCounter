function [layerpos,logPobs] = viterbi(T,nLayerMax,layer0_pos,d,dmax,D,logpd,logb,tiepoints,lastlayerpos,plotlevel)

%% viterbi(T,nLayerMax,layer0_pos,d,dmax,D,logpd,logb,tiepoints,plotlevel):
% delta(t,j,d) represents the maximum likelihood of ending state (j,d) at
% time t, and observing the partial observation sequence obs(1:t):
% delta(t,j,d) = ...
%     max_(s(1):s(t-d)) {P(s(1):s(t-d), s([t-d+1,t])=j, obs(1:t) | theta)}
% where theta signifies all model parameters.
% Also calculated is:
% delta_bar(t,j) = ...
%     max_(s(1):s(t-1)) {P(s(1):s(t-1), s(t])=j, obs(1:t) | theta)}

% The variable psi(t,j) is a back-pointer, showing which state we were
% most likely to arrive from, and used for the subsequent backtracing
% procedure. It records the most likely duration index of state j, assuming
% this state to end at time t.

% Variables:
% T: length of observational data
% nLayerMax: upper-bound estimate of years in observation sequence (fixed in case of tiepoints)
% layer0_pos: probability for ending of previous layer [pixels before obs]
% d: possible durations
% dmax: maximum duration
% D: number of possible durations
% logpd: corresponding log-probabilities
% logb: log-likelihood of data segments to span exactly one year
% tiepoints: depth and age of possible tiepoints

% Copyright (C) 2015  Mai Winstrup
% 2014-10-10 16:15: Clean-up

%% Initialization:
nLayerMax = nLayerMax+1; % Including data for layer 0
logdelta = -ones(T+2*dmax,nLayerMax,D)*Inf; % log(0)=-inf
logdelta_bar = -ones(T+2*dmax,nLayerMax)*Inf;
psi = nan(T+2*dmax,nLayerMax);

% Value for layer 0:
logp0 = log(layer0_pos);
for i = 1:dmax
    logdelta(i,1,:)=logp0(i)+logpd;
end
logdelta_bar(1:dmax,1)=logp0;
    
%% Recursive evaluation of delta:
for t = 1:T
    for j = 2:nLayerMax
        logdelta(t+dmax,j,:)=logdelta_bar(t+dmax-d,j-1)+logb(t,:)'+logpd';
        [logdelta_bar(t+dmax,j), psi(t+dmax,j)]=max(logdelta(t+dmax,j,:));
    end
end
% Note, that the definition of psi here is slightly different to the
% original version. It is transformed such as to reduce the number of
% entries in psi. Its value records the most likely layer duration index
% for the present layer, assuming it to end at t.

% Extending delta for t>T:
for t = T+1:T+dmax
    for j = 2:nLayerMax
        bstart = t-d+1;
        mask = bstart<=T;
        logb2 = zeros(D,1); % For bstart>T+dmax: P(obs)=1 (and log(1)=0)
        logb2(mask)=logb(t,mask);
        
        logdelta(t+dmax,j,:)=logdelta_bar(t+dmax-d,j-1)+logb2+logpd';
        [logdelta_bar(t+dmax,j), psi(t+dmax,j)]=max(logdelta(t+dmax,j,:));
    end
end

% Removing the initialization part of delta and psi, such that these now
% are indexed by the usual t and j:
logdelta = logdelta(dmax+1:end,2:end,:);
logdelta_bar = logdelta_bar(dmax+1:end,2:end);
psi = psi(dmax+1:end,2:end);
% Data for layer 0 is also removed.

if plotlevel>2
    figure;
    plot(logdelta_bar)
    title('log(\delta(t,j))','fontweight','bold')
    xlabel('Pixel')
    ylabel('Prob')
    
    figure;
    pcolor(d(psi'))
    shading flat
    title('d(\psi)','fontweight','bold')
    xlabel('Pixel')
    ylabel('Layer')
    colorbar
end

%% Final maximum likelihood state:
maxprb = nan(1,D);
maxi = nan(1,D);

if isempty(tiepoints)
    maxj = nan(1,D);
    for id = 1:D
        x=logdelta(T:T+d(id)-1,:,id);
        [maxprb(id), index] = max(x(:));
        maxj(id) = ceil(index/d(id)); % ML j-value
        maxi(id) = index-(maxj(id)-1)*d(id); % ML i-value (indexed relative to T)
    end
    maxt = T-1+maxi; % Converting from index i to t
    [logPobs, maxid] = max(maxprb);
    % We have here actually calculated logPobs for the Viterbi algorithm.

    % Most likely layer number, ending time and duration of last layer covered
    % by observation sequence:
    jfinal = maxj(maxid);
    layerpos(jfinal,1) = maxt(maxid);
    dur = d(maxid);
    
elseif ~isempty(lastlayerpos)
    % Final layer is known to be:
    jfinal = nLayerMax-1; 
    % Vi kender placering af det sidste lag fra FB:
    layerpos(jfinal,1) = lastlayerpos;
    
    % Most likely duration for this location:
    x = logdelta(lastlayerpos,jfinal,:);
    [logPobs, maxid] = max(x);
    dur = d(maxid);
    
else
    % With tiepoints:
    % Final layer is known to be:
    jfinal = nLayerMax-1; 
    
    % Most likely duration of layer jfinal:
    for id = 1:D
        % Probability that layer jfinal ends at pixels T:T+dmax having had
        % duration d:
        x = logdelta(T:T+d(id)-1,jfinal,id);
        % Most likely ending position (index) for this value of d:
        [maxprb(id), maxi(id)] = max(x);
    end

    % Most likely duration:
    [logPobs,maxid] = max(maxprb);
    dur = d(maxid);
    % Most likely ending time: Convert from index to pixels:
    maxt = T-1+maxi(maxid);
    layerpos(jfinal,1) = maxt;   

end

%% Backtracking to find the maximum likelihood sequence:
for j = jfinal:-1:2
    % Previous layer ends at time:
    layerpos(j-1,1) = layerpos(j)-dur;
    % Having a most likely duration of:
    dur = d(psi(layerpos(j-1,1),j-1));
end
% psi(t,j) records the most likely duration index of state j, assuming this
% state to end at time t.

% Removing the last layer boundary if it occurs after end of the data
% sequence:
layerpos(layerpos>T)=[];

%% All layer boundaries are moved to be located midway between data points:
layerpos = layerpos+0.5;