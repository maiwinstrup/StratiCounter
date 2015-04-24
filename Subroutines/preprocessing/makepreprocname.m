function preprocname = makepreprocname(preprocsteps,dx)

%% preprocname = makepreprocname(preprocsteps,dx)
% This function creates a file/folder name corresponding to a specific 
% data preprocessing and resolution.
% Copyright (C) 2015  Mai Winstrup

%% Name for preprocessing steps:
preprocname = [];
for i = 1:size(preprocsteps,1)
    % Numeric values (distances etc.) used in preprocessing:
    if size(preprocsteps,2)==1
        specs = [];
    else
        % Distance value:
        dists = num2str(preprocsteps{i,2});
        % Additional specifications:
        values = [];
        if size(preprocsteps(i,:),2)>2
            numval = preprocsteps{i,3};
            for k=1:length(numval); 
                values = [values ',' num2str(numval(k))]; 
            end
        end
        if isempty(dists)
            specs = values(2:end);
        elseif isempty(values)
            specs = dists;
        else
            specs = [dists ',' values];
        end
    end
    % Append to filename:
    preprocname = [preprocname preprocsteps{i} specs '_'];
end
preprocname = preprocname(1:end-1);

%% Add resolution of data file:
% Number of significant digits:
if dx>=10^-2; 
    res = [num2str(dx*100) 'cm'];
else res = [num2str(dx*10^3) 'mm'];
end
if isempty(preprocname); preprocname = ['none_' res];
else preprocname = [preprocname '_' res];
end