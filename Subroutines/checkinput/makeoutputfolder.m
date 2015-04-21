function [outputdir,outputdir0,runID] = makeoutputfolder(Model,Runtype)

%% [outputdir,outputdir0,runID] = makeoutputfolder(Model,Runtype)
% Generate outputfolders for results. Outputdir0 is directory corresponding
% to current model settings, while outputdir also includes runID. 
% Copyright (C) 2015  Mai Winstrup

%% Find/generate output folder:
% Running in development mode?
if strcmp(Runtype.develop,'yes'); outputdir0 = './Output/develop';
else outputdir0 = './Output';
end

if strcmp(Model.icecore, 'SyntheticData')
    outputdir0 = [outputdir0 '/SyntheticData/' Model.info];
else 
    outputdir0 = [outputdir0 '/' Model.icecore '/Timescale/' ...
        num2str(Model.dstart) '-' num2str(Model.dend) 'm'];
end

%% Retrieve runID:
if exist([outputdir0 '/runID.mat'],'file')
    load([outputdir0 '/runID.mat']);
    runID = runID+1;
else
    runID = 1;
end

%% Output directory for current run:
outputdir = [outputdir0 '/#' num2str(runID)];

%% Make folder if doesn't exist, and remove its content if it does:
if ~exist(outputdir,'dir')
    mkdir(outputdir); 
else
    % Remove content in folder from old run:
    disp(['Deleting content in ' outputdir])
    delete([outputdir '/*.*'])
end