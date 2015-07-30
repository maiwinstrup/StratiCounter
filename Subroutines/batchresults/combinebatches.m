function [Layerpos,LayerProbDist,centralEst,timescale,timescale1yr,markerProb,...
    markerConf,lambda] = combinebatches(Result,manualcounts,Model)

%% [Layerpos,LayerProbDist,centralEst,timescale,timescale_1yr,markerProb,...
%       markerConfInt,lambda] = combinebatches(Result,manualcounts,Model)
% Combining results from the individual batches, converting to timescales
% using the selected age unit (ageUnitOut), computing results from number
% of layers between marker horizons, computing mean layer thicknesses. 

% Output: 
% Layerpos: Layer positions directly from batch results
%           Layerpos.fb, Layerpos.fb_issues: Forward-Backward
%           Layerpos.viterbi: Viterbi results, Layerpos.combined = Improved 
%           location of layer boundaries; found using viterbi algorithm 
%           when constrained by the Forward-Backward results. 
% LayerProbDist: Layer probability distribution in "ageUnitOut" for all data 
%           LayerProbDist.d: Depth, LayerProbDist.mode, LayerProbDist.mean, 
%           LayerProbDist.median, LayerProbDist.prctile (all from FB results)
% centralEst: [depth, mode, mean, median] for all depths, in "ageUnitOut"
%           (derived from FB results).
% timescale: [depth, mode, percentiles] for all depths, in "ageUnitOut"
%           (derived from FB results).
% timescale1yr: [depth_layerboundary, age, percentiles] in "ageUnitOut"
%           Layer boundaries are those found using the viterbi algorithm 
%           and constrained by the Forward-Backward results. 
%           Percentiles of the confidence interval is given at midpoints of 
%           layers in order to eliminate uncertainties related to the
%           location of the layer boundary. 
% markerProb{iMarkerSet}: Probability distributions for marker horizons
%           markerProb{iMarkerSet}(imarker).ndist = [#layer, prob]
%           markerProb{iMarkerSet}(imarker).d = depth (end of interval)
% markerConf{iMarkerSet}: Mode, probability of mode, and confidence 
%           intervals for number of layers between successive marker horizons
%           [dstart, dend, #layer_ML, prob(#layer_ML), confidence intervals]
% lambda{idx}: Layer thicknesses at various equidistant intervals 
%           [dstart, dend, lambda_ML, confidence interval for lambda]

% Copyright (C) 2015  Mai Winstrup

%% Number of batches:
nBatch = length(Result);

%% Layer boundaries directly from batch output:
% For Forward-Backward and the Viterbi algorithms: 
Layerpos.fb = [];
Layerpos.fb_issues = [];
Layerpos.final = [];
for iBatch = 1:nBatch
    Layerpos.fb = [Layerpos.fb; Result(iBatch).Layerpos.fb];
    Layerpos.fb_issues = [Layerpos.fb_issues; Result(iBatch).Layerpos.fb_issues];
    Layerpos.final = [Layerpos.final; Result(iBatch).Layerpos.final];
end

% Give warning if the algorithm encountered problems converting from the
% full timescale probability distributions to layerboundaries:
if ~isempty(Layerpos.fb_issues)
    disp(['Warning: Issue occurred during placement of FB layer'...
        ' boundaries at the following depths:'])
    for i = 1:size(Layerpos.fb_issues,1)
        disp([num2str(Layerpos.fb_issues(i,1)) 'm: ' ...
            num2str(Layerpos.fb_issues(i,2)) ' layer(s)'])
    end
    disp(['Check that the layering at these locations have been resolved'...
        ' in a satisfying manner'])
end

%% Full results from Forward-Backward algorithm:
% Probability distributions are summarized using their mode, median and 
% mean, and confidence intervals are provided using percentiles of the 
% layer distributions for each data point. 
LayerProbDist = struct('d',[],'mode',[],'mean',[],'median',[],'prctile',[]);
for iBatch = 1:nBatch
    LayerProbDist.d = [LayerProbDist.d; Result(iBatch).LayerProbDist.d(1:end-1)];
    LayerProbDist.mode = [LayerProbDist.mode; Result(iBatch).LayerProbDist.mode(1:end-1)];
    LayerProbDist.mean = [LayerProbDist.mean; Result(iBatch).LayerProbDist.mean(1:end-1)];
    LayerProbDist.median = [LayerProbDist.median; Result(iBatch).LayerProbDist.median(1:end-1)];
    LayerProbDist.prctile = [LayerProbDist.prctile; Result(iBatch).LayerProbDist.prctile(1:end-1,:)];
end
% And the last data point:
LayerProbDist.d = [LayerProbDist.d; Result(iBatch).LayerProbDist.d(end)];
LayerProbDist.mode = [LayerProbDist.mode; Result(iBatch).LayerProbDist.mode(end)];
LayerProbDist.mean = [LayerProbDist.mean; Result(iBatch).LayerProbDist.mean(end)];
LayerProbDist.median = [LayerProbDist.median; Result(iBatch).LayerProbDist.median(end)];
LayerProbDist.prctile = [LayerProbDist.prctile; Result(iBatch).LayerProbDist.prctile(end,:)];

%% Convert to ages instead of layer number:
% Estimate/find age of the first autocounted annual layer boundary 
% (according to the manual counts):
t0 = initialage(manualcounts,[LayerProbDist.d(:),LayerProbDist.mode(:)],Model);
% While t0 can often be inferred directly, only an estimate can be provided 
% if manual counts do not cover the first layer boundary. 
        
% Pixel 1 is part of "layer 0", and the corresponding age here is t0 
% (measured in appropriate age unit). 
switch Model.ageUnitOut
    case 'AD'
        LayerProbDist.mode = t0-LayerProbDist.mode;
        LayerProbDist.mean = t0-LayerProbDist.mean;
        LayerProbDist.median = t0-LayerProbDist.median;
        LayerProbDist.prctile = t0-LayerProbDist.prctile;
    otherwise
        LayerProbDist.mode = t0+LayerProbDist.mode; 
        LayerProbDist.mean = t0+LayerProbDist.mean;
        LayerProbDist.median = t0+LayerProbDist.median;
        LayerProbDist.prctile = t0+LayerProbDist.prctile;   
end

%% Combining results into various arrays:
% Central estimates: Mode, mean, median
centralEst = [LayerProbDist.d, LayerProbDist.mode, LayerProbDist.mean, LayerProbDist.median];

% Combined timescale output: 
% Mode and confidence intervals (percentiles of distribution)
timescale = [LayerProbDist.d, LayerProbDist.mode, LayerProbDist.prctile];

%% Produce timescale array in 1-yr resolution: 
% The layer boundaries are found using a combined approach: The viterbi 
% algorithm is employed to find the most likely layer boundaries, when
% constrained by the total most likely number of layers (for each batch) as 
% derived from the Forward-Backward results. We then know that each of the
% boundaries should be counted as exactly one year. 
% Confidence intervals are based on percentiles of the layer number 
% probability distributions, and taken midway between layer boundaries in 
% order to eliminate uncertainties due to the exact boundary locations. 

% Depth of boundaries, incl. initial layer position:
timescale1yr(:,1) = [LayerProbDist.d(1); Layerpos.final];

% Corresponding ages, incl. that of initial layer position:
switch Model.ageUnitOut
    case 'AD'
        timescale1yr(:,2) = t0-(-1:length(Layerpos.final)-1);
    otherwise
        timescale1yr(:,2) = t0+(0:length(Layerpos.final));
end

% Confindence interval estimates are taken at midpoints of layers:
% In this way, we remove uncertainties related to the exact placement of 
% the layer boundaries.
layerpos_px = interp1(LayerProbDist.d,1:length(LayerProbDist.d),...
    [LayerProbDist.d(1);Layerpos.final-0.5*Model.dx]); % Last pixel in layers
switch Model.ageUnitOut
    case 'AD'
        layercenter_px = round(layerpos_px-0.5*diff([1; layerpos_px]));
    otherwise 
        layercenter_px = round(layerpos_px+0.5*diff([layerpos_px; length(LayerProbDist.d)]));
end
% Confidence interval(s):
timescale1yr(:,3:2+size(LayerProbDist.prctile,2)) = ...
    LayerProbDist.prctile(layercenter_px,:);
% Correct for initial layer boundary:
timescale1yr(1,3:2+size(LayerProbDist.prctile,2)) = timescale1yr(1,2);

%% Layer number distributions between marker horizons:
% Several sets of marker horizons may exist, the probabilities of these are
% calculated one set at a time.
markerProb = cell(1,length(Model.dMarker));
zerolimit = 10^-5;
for iMarkerSet = 1:length(Model.dMarker)
    k=0;
    for iBatch = 1:nBatch
        % Number of marker horizons within batch (end of sections):
        nMarker = length(Result(iBatch).Marker(iMarkerSet).d);
        for j = 1:nMarker
            % End depth:
            markerProb{iMarkerSet}(k+j).d = Result(iBatch).Marker(iMarkerSet).d(j);
            % Distribution:
            markerProb{iMarkerSet}(k+j).ndist = compactprobdist(...
                [Result(iBatch).Marker(iMarkerSet).ndist(:,1), ...
                Result(iBatch).Marker(iMarkerSet).ndist(:,j+1)],zerolimit);
        end
        k = k+nMarker;
    end
end

%% Maximum likelihood and confidence intervals of layer number probability
% distributions between marker horizons: 
% Also the section from the beginning to the first marker horizon is 
% included, as well as the section from the last marker horizon within data 
% section to the last data point.
markerConf = cell(1,length(Model.dMarker));

for iMarkerSet = 1:length(Model.dMarker)
    % Remove marker horizons outside interval:
    mask = Model.dMarker{iMarkerSet}>=LayerProbDist.d(1) & ...
        Model.dMarker{iMarkerSet}<=LayerProbDist.d(end);
    markerhorizons = Model.dMarker{iMarkerSet}(mask); % Includes start of 
    % data section.
    % Number of horizons within interval:
    nMarker = length(markerhorizons); 
    if nMarker <= 1; continue; end
    
    % Corresponding depth intervals: 
    dstart = markerhorizons(1:nMarker-1);
    dend = markerhorizons(2:nMarker);
    
    % Initialize:
    probML = nan(nMarker-1,1);
    nML = nan(nMarker-1,1);
    confInt = nan(nMarker-1,length(Model.prctile));
    
    for iMarker = 1:nMarker-1
        % Maximum likelihood layer number, and its associated probability:
        [probML(iMarker),indexML] = ...
            max(markerProb{iMarkerSet}(iMarker).ndist(:,2));
        nML(iMarker) = markerProb{iMarkerSet}(iMarker).ndist(indexML,1);
        % Confidence intervals:
        confInt(iMarker,:) = ...
            prctileofprobdist(markerProb{iMarkerSet}(iMarker).ndist(:,1),...
            markerProb{iMarkerSet}(iMarker).ndist(:,2),Model.prctile)'; 
    end
    
    % Combine to array:
    markerConf{iMarkerSet} = [dstart(:), dend(:), nML, probML, confInt];
end

%% Mean layer thicknesses in equidistant intervals:
lambda = cell(1,length(Model.dxLambda));

for idx = 1:length(Model.dxLambda)    
    % Combine results from batches:
    nMode = []; prc = []; dstart = LayerProbDist.d(1);        
    for iBatch = 1:nBatch
        if isempty(Result(iBatch).Lambda(idx).ndist)
            continue
        end
        
        % Layer number probability distributions:
        n = Result(iBatch).Lambda(idx).ndist(:,1);
        p = Result(iBatch).Lambda(idx).ndist(:,2:end);
        
        % Calculating maximum likelihood layer numbers in intervals
        % This is mode of the probability distributions:
        [~,index]=max(p,[],1); 
        nMode = [nMode; n(index(:))];
            
        % Calculating percentiles:
        prc = [prc; prctileofprobdist(n,p,Model.prctile)'];
        
        % Start depth of intervals:
        dstart = [dstart; Result(iBatch).Lambda(idx).d(:)];
    end
    
    % Construct lambda array for this value of dxlambda:
    if ~isempty(nMode)
        % Interval:
        lambda{idx}(:,1) = dstart(1:end-1);
        lambda{idx}(:,2) = dstart(2:end);
    
        % Convert from layer numbers to layer thicknesses and their 
        % associated confidence interval:
        % Size of the individual intervals:
        L = diff(dstart);
        % Most likely mean layer thicknesses:
        lambda{idx}(:,3) = L./nMode;
        % Confidence intervals for mean layer thickness (from smallest to 
        % largest value):
        nPrc = length(Model.prctile);
        lambda{idx}(:,4:3+nPrc) = repmat(L,1,nPrc)./fliplr(prc); 
    end
end