function [preprocsteps,preprocdist] = setpreprocdist(preprocsteps0,meanLambda)

%% [preprocsteps,preprocdist] = setpreprocdist(preprocsteps0,meanLambda)
% The preprocessing distances in the array preprocsteps0 are assumed as 
% given in units of mean layer thicknesses (meanLambda). These values are 
% here converted to units of meters. An adjusted preprocessing array is
% provided (preprocsteps) along with the numeric values of the employed 
% preprocessing distances (preprocdist).
% Copyright (C) 2015  Mai Winstrup

%% Make preprocessing array with adjusted distances:
nSpecies = size(preprocsteps0,1);
preprocsteps = cell(nSpecies,1);
preprocdist = cell(nSpecies,1);

% Processing types:
preprocsteps(:,1) = preprocsteps0(:,1);

% Processing distances, rounded to two significant digits:
for j = 1:nSpecies
    if ~isempty(preprocsteps0{j})
        preprocdist{j,1} = roundsignificant(...
            cell2mat(preprocsteps0{j}(:,2))*meanLambda,2); %[m]
        for k = 1:length(preprocdist{j}(:,1))
            preprocsteps{j,1}(k,2)={preprocdist{j}(k)};
        end
    end
end

% Additional processing specs (if any):
for j = 1:nSpecies
    if size(preprocsteps0{j,1},2)>2
        preprocsteps{j,1}(1,3)= preprocsteps0{j,1}(1,3);
    end
end