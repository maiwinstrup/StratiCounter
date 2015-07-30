function [Result, Layer0_new, batchStart_new] = ...
    resultsforbatch(depth,FBprob,Layerpos,Layer0,d,pd,logb,meanLambda,...
    batchStart,dDxLambda,iIter,Model,plotlevel)

%% [Result, Layer0_new, batchStart_new] = resultsforbatch(depth,FBprob,...
%    Layerpos,Layer0,d,pd,logb,meanLambda,batchStart,dDxLambda,iIter,...
%    Model,plotlevel)
% Calculating the results (layer probability distribution and most likely 
% layer boundaries) for the current batch, and setting initial conditions 
% for the next. Also calculated are probability distributions between 
% marker horizons, and for regular intervals (later used for calculating 
% mean layer thicknesses within these). 

% Output:
% batchStart_new: Start pixel for next batch
% Layer0_new.pos: Probability of where the next "layer0" is starting 
% Layer0_new.no: Layer number distribution at tau 
% Layer0_new.noDx, Layer0_new.noMarker: Same, but for layer thickness 
% intervals and for marker horisons.  
% Result.LayerDist.d: depth scale
% Result.LayerDist.mode: mode of timescale
% Result.LayerDist.mean: mean of timescale
% Result.LayerDist.quantile: quantiles of timescale
% Result.Layerpos.fb: layer positions according to FB algorithm [m]
% Result.Layerpos.final: a "best" set of layer boundaries (using FB and 
% Viterbi algorithms combined)
% Result.Lambda.ndist: Probability distribution of layers in section 
% Result.Lambda.d: ending depth of corresponding lambda sections
% Result.Marker.ndist: Probability distribution of layer numbers in section
% Result.Marker.d: End depth of section for marker horizons.
% Result.nIter: Final iteration number for batch

% Copyright (C) 2015  Mai Winstrup
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.

%% Select ending pixel for batch: 
% We only consider the first part of data sequence (i.e. to and including 
% pixel tau), as the annual layering in the last part can be better 
% estimated when also using the subsequent data.
% The pixel tau is selected as the first pixel after a very likely layer 
% boundary, for which we are sure to be out of the previous layer. 
% Also calculating the ending location of the layer previous to the one in 
% pixel tau. 
batchLength = size(FBprob.gamma,1);
[tau, postau] = findtau(FBprob,meanLambda,d,batchLength,Model,plotlevel); 
% Note: tau is included in both batches

% Use for initialization of next batch:
Layer0_new.pos = postau;
% Next batch starts at:
batchStart_new = batchStart+tau-1; %[pixel from start of data series]

%% Resulting layer number distribution along batch:
% This is calculated up to and including pixel tau.
[ntau, ntauTotal, LayerProbDist] = ...
    batchlayerprobabilities(depth,FBprob,tau,Layer0,Model.prctile,plotlevel); 

% Save: Results for current batch, and initial conditions for next
Layer0_new.no = ntauTotal;
Result.LayerProbDist = LayerProbDist;

%% The most likely layer boundaries:
% Converting to depth, and computing an optimal set of layer boundaries by 
% combining results from the Forward-Backward algorithm with the Viterbi
% approach.
% OBS: Layerpos.final has not been checked!
[LayerposDepth, ~] = batchlayerpos(Layerpos,depth,...
    tau,Layer0,postau,ntauTotal,d,pd,logb,Model.dx,plotlevel);
Result.Layerpos = LayerposDepth;

%% Mean layer thickness within sections in batch:
% Layer number probability distributions for predetermined sections; these 
% are later to be used for calculating mean layer thickness:
for ix = 1:length(Model.dxLambda)
    % Layer number distributions within sections:
    [probdistSections,Layer0_new.noDx{ix},dSectionBounds] = ...
        sectionlayerprobs(depth,FBprob,dDxLambda{ix},Layer0.noDx{ix},...
        tau,ntau,d,pd,logb,plotlevel);
   
    % Save: 
    Result.Lambda(ix).d=[];
    Result.Lambda(ix).ndist=[];
    if ~isempty(probdistSections)
        Result.Lambda(ix).d = dSectionBounds; 
        Result.Lambda(ix).ndist = probdistSections;
    end
end

%% Marker horizons: 
% Layer number probability distributions for sections between marker
% horizons:
for ix = 1:length(Model.dMarker)
    % Layer number distributions within sections:
    [probdistMarker,Layer0_new.noMarker{ix},dMarkerBounds] = ...
        sectionlayerprobs(depth,FBprob,Model.dMarker{ix},...
        Layer0.noMarker{ix},tau,ntau,d,pd,logb,plotlevel);
   
    % Save: 
    Result.Marker(ix).d = dMarkerBounds; 
    Result.Marker(ix).ndist = probdistMarker;
end

%% Number of iterations performed:
% Total iteration number: 
Result.nIter = iIter;