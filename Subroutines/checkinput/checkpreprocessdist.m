function checkpreprocessdist(Model,manualcounts)

%% checkpreprocessdist(Model,manualcounts)
% Check that the values of model.dx and preprocessing distance are chosen 
% reasonably compared to the layer thicknesses in interval.
% The following generally seems to work well: 
% dx < 10*lambda, dist > 2*lambda

% Mai Winstrup
% 2014-06-17 21:09: Initial version
% 2014-06-19 17:00: Only performed when given manual counts
% 2014-10-01 12:52: Floating distances accounted for separately

%% Only performed if manual counts are provided:
if isempty(manualcounts); return; end

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
        disp(['OBS: Average number of data points per layer is small (' num2str(Naverage) ')'])
    end
    if Naverage > 20
        disp(['OBS: Average number of data points per layer is large (' num2str(Naverage) ')'])
    end
    if Nmin < 5
        disp(['OBS: Minimum number of data points per layer: ' num2str(Nmin)])
    end
end

%% Preprocessing distance (minimum value) relative to the layer thickness:
minDist = inf;
minDistFloat = inf;
for j = 1:Model.nSpecies
    preproc_type = Model.preprocess{j}(:,1);
    stringpos = strfind(preproc_type,'_float');
    dists = Model.preprocess{j}(:,2);
    for k = 1:length(dists)
        dist_here = cell2mat(dists(k));
        if isempty(stringpos{k}) % Not floating dist
            minDist = min(minDist,dist_here);
            % Select minimum value encountered
        else
            minDistFloat = min(minDistFloat,dist_here);
        end
    end
end

if minDistFloat < 2
    disp('OBS: Floating distance is less than 2*lambda')
end

if isfinite(minDist)
    % Average number of layers per this distance:
    Maverage = minDist/meanLambda;
    % Minimum number of layers per this distance:
    Mmin = minDist/maxLambda;
    disp(['Average number of layers per (smallest) preprocessing distance: ' num2str(Maverage)])
    if Maverage < 2
        disp(['OBS: Average number of layers per preprocessing distance is small (' num2str(Maverage) ')'])
    end
    if Mmin < 1
        disp(['OBS: Minimum number of layers per preprocessing distance is small: ' num2str(Mmin)])
    end
end