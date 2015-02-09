function [batchStart,Layer0,Template,Prior,Layerpar,logPobs,relweight,Result]=...
    initializematrices(nBatch,Model)

%% [batchStart,Layer0,Template,Prior,Layerpar,logPobs,relweight,Result]=...
%    initializematrices(nBatch,Model)
% Initialize empty matrices for results from nBatches. 

% Copyright (C) 2015  Mai Winstrup
% 2014-10-22 14:37

%% Start pixel for batches:
batchStart = nan(nBatch,1);

%% Probabilities for location and layer number of first layer boundary:
if isempty(Model.dMarker)
    Layer0(1:nBatch) = struct('pos',[],'no',[],'noDx',[]);
else
    Layer0(1:nBatch) = struct('pos',[],'no',[],'noDx',[],'noMarker',[]);
end

%% Initialize layer templates and layer parameters:
switch Model.type
    case 'PCA'
        Template(1:Model.nSpecies,1:nBatch,1:Model.nTemplateBatch+1) = ...
            struct('mean',[],'dmean',[],'d2mean',[],'traj',[],'dtraj',[],'d2traj',[]);
    case 'FFTcomp'
        Template(1:Model.nSpecies,1:nBatch,1:Model.nTemplateBatch+1) = ...
            struct('dc',[],'phase',[],'amplitude',[]);
end

%% Prior on parameter estimates:
Prior(1:nBatch,1:Model.nTemplateBatch) = ...
    struct('m',[],'v',[],'sigma',[],'u',[],'invU',[],'cov',[],'nvar',[]);

%% Layer parameters:
Layerpar(1:nBatch,1:Model.nTemplateBatch,1:Model.nIter)=...
    struct('my',nan,'sigma',nan,'par',nan(Model.nSpecies,Model.order),...
    'cov',nan(Model.nSpecies*Model.order,Model.nSpecies*Model.order),...
    'nvar',nan(Model.nSpecies,1)); 

%% Initialize matrices for performance evaluation: 
% Evolution of log(Pobs) with iterations:
logPobs = nan(nBatch,Model.nTemplateBatch,Model.nIter+1,2);

% Relative weighting of layer shape vs. layer thickness:
relweight = nan(nBatch,Model.nTemplateBatch,Model.nIter);     

%% Initialize matrices containing batch results: 
% Probability results:
LayerDist = struct('d',[],'mode',[],'mean',[],'median',[],'prctile',[]);
% Layer positions:
Layerpos = struct('fb',[],'fb_issues',[],'viterbi',[],'combined',[]);
% Mean layer thickness within intervals:
Lambda = struct('d',[],'ndist',[]);
Lambda = repmat(Lambda,length(Model.dxLambda),1);
% Probability distributions of layering between marker horizons:
Marker = struct('d',[],'ndist',[]);
Marker = repmat(Marker,length(Model.dMarker),1);
% Gathering all in the array "Result": 
if isempty(Model.dMarker)
    Result(1:nBatch) = struct('LayerDist',LayerDist,'Layerpos',Layerpos,...
        'Lambda',Lambda,'nIter',[]);
else
    Result(1:nBatch) = struct('LayerDist',LayerDist,'Layerpos',Layerpos,...
        'Lambda',Lambda,'Marker',Marker,'nIter',[]);
end