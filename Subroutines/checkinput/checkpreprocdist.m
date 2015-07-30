function checkpreprocdist(Model,manualcounts)

%% checkpreprocdist(Model,manualcounts)
% Check that the values of model.dx and preprocessing distances are chosen 
% reasonably compared to the layer thicknesses in interval.
% The following generally seems to work well: 
% dx < 10*lambda, dist > 2*lambda

% Copyright (C) 2015  Mai Winstrup

%% Only performed if manual counts are provided:
if isempty(manualcounts); 
    mindist = [];
    for j = 1:Model.nSpecies
        mindist = min([mindist; cell2mat(Model.preprocsteps{j,1}(:,2))]);
    end
    if ~isempty(mindist)
        warning(['Check that preprocessing distances (min value: ' ...
            num2str(mindist) ') are much larger than derived layer thicknesses'])
    end
    return
end

%% Range of layer thicknesses for the interval [m]: 
meanLambda = mean(diff(manualcounts(:,1)));
minLambda = min(diff(manualcounts(:,1)));
maxLambda = max(diff(manualcounts(:,1)));

%% Number of data points per layer:
if ~isempty(Model.dx)
    Naverage = round(meanLambda/Model.dx);
    Nmin = round(minLambda/Model.dx);
    disp(['Average number of data points per layer: ' num2str(Naverage)])
    if Naverage < 10
        warning(['OBS: Average number of data points per layer is small (' num2str(Naverage) ')'])
    end
    if Naverage > 20
        warning(['OBS: Average number of data points per layer is large (' num2str(Naverage) ')'])
    end
    if Nmin < 5
        warning(['OBS: Minimum number of data points per layer: ' num2str(Nmin)])
    end
end

%% Preprocessing distance (minimum value) relative to the layer thickness:
% Fixed preprocessing distances:
preprocdist{Model.nSpecies,1}=[];
for j = 1:Model.nSpecies
    if ~isempty(Model.preprocsteps{j,1}) && size(Model.preprocsteps{j,1},2)>1
        preprocdist{j,1} = cell2mat(Model.preprocsteps{j,1}(:,2)); %[m]
    end
end
minDist = min(cell2mat(preprocdist));

if isfinite(minDist)
    % Average number of layers per this distance:
    Maverage = minDist/meanLambda;
    % Minimum number of layers per this distance:
    Mmin = minDist/maxLambda;
    disp(['Average number of layers per (smallest) preprocessing distance: ' num2str(Maverage)])
    if Maverage < 2
        warning(['OBS: Average number of layers per preprocessing distance is small (' num2str(Maverage) ')'])
    end
    if Mmin < 1
        warning(['OBS: Minimum number of layers per preprocessing distance is small: ' num2str(Mmin)])
    end
end

% Adjustable preprocessing distances: 
preprocdist_adj{Model.nSpecies,1}=[];
for j = 1:Model.nSpecies
    if ~isempty(Model.preprocsteps{j,2}) && size(Model.preprocsteps{j,2},2)>1
        preprocdist_adj{j,1} = cell2mat(Model.preprocsteps{j,2}(:,2)); %[fractions of lambda]
    end
end
minDist_adj = min(cell2mat(preprocdist_adj));

if minDist_adj < 2
     warning('OBS: Adjustable preprocessing distance(s) is small (less than 2*lambda)')
end