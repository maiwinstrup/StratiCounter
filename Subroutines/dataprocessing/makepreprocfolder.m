function outputdir = makepreprocfolder(preprocess,dx,icecore,species,Runtype)

%% outputdir = makepreprocfolder(preprocess,dx,icecore,species,Runtype)
% This function creates a folder/foldername corresponding to a specific 
% data processing.

% Copyright (C) 2015  Mai Winstrup
% 2014-04-21 17:01: Filename updated
% 2014-07-15 14:11: New structure of preproc
% 2014-08-15 16:24: Small adjustments
% 2014-08-16 13:37: New name of function
% 2014-08-21 21:12: Function updated, and converted to make corresponding 
%                   folder/foldername

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

%% Output folder for processed data:
outputdir = [outputdir '/' icecore '/ProcessedData/' species '/' preprocname];

%% Make directory:
if ~exist(outputdir,'dir')
    mkdir(outputdir)
end