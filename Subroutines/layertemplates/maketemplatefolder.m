function outputdir = maketemplatefolder(preprocess,dx,icecore,species,type,normalizelayer,Runtype)

%% outputdir = maketemplatefolder(preprocess,dx,icecore,species,type,normalizelayer,Runtype)
% This function creates a folder/foldername corresponding to a specific 
% data processing and the employed layer model.

% Copyright (C) 2015  Mai Winstrup
% 2014-08-21 21:19: Initial verison
% 2014-08-22 15:00: Minor changes

%% Running in development mode?
if strcmp(Runtype.develop,'yes'); outputdir = './Output/develop';
else outputdir = './Output';
end

%% Including the (ordered) preprocessing steps:
preprocname = [];
for i = 1:size(preprocess,1)
    % Distance(s) used for this processing:
    valnum = preprocess{i,2};

    % Convert to string:
    values = [];
    for k=1:length(valnum)
        values = [num2str(valnum) '-'];
    end
    values = values(1:end-1);
    
    % Append to filename:
    preprocname = [preprocname preprocess{i} values '_'];
end
preprocname = preprocname(1:end-1);

%% Add resolution of data file:
% Number of significant digits:
if dx*100 >= 1; 
    res = [num2str(dx*100) 'cm'];
else res = [num2str(dx*10^3) 'mm'];
end
preprocname = [preprocname '_' res];

%% Output folder for manual layer templates and parameters:
outputdir = [outputdir '/' icecore '/LayerCharacteristics/' species '/' preprocname '/' type];
if ~strcmp(normalizelayer,'none')
    outputdir = [outputdir '_' normalizelayer];
end

%% Make directory:
if ~exist(outputdir,'dir')
    mkdir(outputdir)
end