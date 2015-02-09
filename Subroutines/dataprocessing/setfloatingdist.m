function preprocess_new = setfloatingdist(preprocess,meanLambda)

%% preprocess_new = setfloatingdist(preprocess,meanLambda)
% Convert floating preprocessing distances to their appropriate values, as 
% determined from the mean layer thickness. 

% Copyright (C) 2015  Mai Winstrup
% 2014-05-13 01:37: Second output changed from index to floating_dist
%                   Input changed: "interval", "counts" -> "meanlambda"
% 2014-05-18 14:20: Outputname changed to floatdist
% 2014-06-13 13:57: meanlambda -> meanLambda
% 2014-07-21 22:53: Changed to reflect new version of preprocessing array
% 2014-08-16 12:27: Round to two significant digits

%% Replace string and value in preprocessing array:
preprocess_new = preprocess;

for j = 1:length(preprocess)
    preproc_type = preprocess{j}(:,1);
    stringpos = strfind(preproc_type,'_float');
    
    for istep = 1:length(stringpos)
        if ~isempty(stringpos{istep});
            % Change name (delete "_float")
            preprocess_new{j}(istep,1) = {preproc_type{istep}(1:cell2mat(stringpos(istep))-1)};
            % Change processing distance: 
            newdist = cell2mat(preprocess{j}(istep,2))*meanLambda;
            % Round to two significant digits: 
            preprocess_new{j}(istep,2) = {roundsignificant(newdist,2)};
        end
    end
end