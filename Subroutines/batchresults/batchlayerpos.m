function [LayerposDepth, logPobs_combined] = ...
    batchlayerpos(Layerpos,depth,tau,Layer0,postau,ntau,d,pd,logb,Model,plotlevel)

%% [LayerposDepth, logPobs_combined] = batchlayerpos(Layerpos,depth,...
% tau,Layer0,postau,ntau,d,pd,logb,Model,plotlevel)
% Converting the entries in Layerpos from pixel to depth, and computing an
% optimal set of layer boundaries based on a Forward-Backward constrained 
% version of the viterbi algorithm.  
% 
% Mai Winstrup
% 2014-10-15 11:22: Separat script, added optimal layer boundaries

%% ForwardBackward and/or Viterbi layering within considered interval:
% Boundaries located in pixel 1 are defined to be part of the previous 
% batch.

% Layering according to the Forward-Backward algorithm, converted to
% depths:
mask = Layerpos.fb(:,1)>1 & Layerpos.fb(:,1)<=tau;
Layerpos.fb = Layerpos.fb(mask);
LayerposDepth.fb = depth(Layerpos.fb-0.5)+0.5*Model.dx;

% Locations with layer location issues:
mask = Layerpos.fb_issues(:,1)>1 & Layerpos.fb_issues(:,1)<=tau;
Layerpos.fb_issues = Layerpos.fb_issues(mask,:);
LayerposDepth.fb_issues(:,1) = depth(Layerpos.fb_issues(:,1)-0.5)+0.5*Model.dx;
LayerposDepth.fb_issues(:,2) = Layerpos.fb_issues(:,2);

% According to the Viterbi algorithm:
if isempty(Layerpos.viterbi)
    LayerposDepth.viterbi = [];
else
    mask = Layerpos.viterbi>1 & Layerpos.viterbi<=tau;
    Layerpos.viterbi = Layerpos.viterbi(mask);
    LayerposDepth.viterbi = depth(Layerpos.viterbi-0.5)+0.5*Model.dx;
end

%% An optimal set of layerboundaries: 
% Found by using the viterbi algorithm, constrained by output of the
% Forward-Backward algorithm. 
% Constraints are: Probability distribution of location of first layer 
% (layer0), location of the last layer boundary (fixed), and the derived 
% most likely number of layers inbetween. 

% Layer duration parameters:
dmax = max(d);
D = length(d);

% Constraints:
% Probability distribution for ending of layer previous to batch:
layer0_pos = zeros(1,dmax);
istart = dmax-length(Layer0.pos)+1;
layer0_pos(max(1,istart):end)=Layer0.pos(max(1,-istart+2):end);
% Normalize probabilities:
layer0_pos = layer0_pos/sum(layer0_pos); %[probabilities]

% Location of last layer boundary (end of layer) in batch:
% Selected as mode of the probability distribution postau.
[~, imax] = max(postau);
lastlayerpx = tau-length(postau)+imax-1; %[pixel]
% We know that this layer boundary will be located before tau, just as will 
% the corresponding Forward-Backward layer boundary.

% Most likely number of layers in batch (counted from zero):
[~,imax] = max(ntau(:,2));
%nLayerML = ntau(imax,1);
nLayerML = ntau(imax,1)+1;

% Using tiepoints:
tiepoints = '42'; % det er lige meget hvad der er heri, blot ikke er tom

[layerpos_combined, logPobs_combined] = viterbi(tau,nLayerML,layer0_pos,d,dmax,D,...
    log(pd),logb(1:tau+dmax,:),tiepoints,lastlayerpx,plotlevel); % 2014-04-02 14:52

% [pixel] - senere!% bør tjekke om det ikke bør være plus dx/2 istedet. er det start eller slutpunkter vi har her? 

% Convert to depth:
LayerposDepth.combined = depth(layerpos_combined-0.5)+0.5*Model.dx;

%% Compare resulting layer boundary positions:
if plotlevel > 1
    figure;
    plot(LayerposDepth.fb,'-k')
    hold on
    plot(LayerposDepth.viterbi,'-b')
    plot(LayerposDepth.combined,'-r')
end