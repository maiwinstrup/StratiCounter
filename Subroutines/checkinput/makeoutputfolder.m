function [outputdir,outputdir0,runID] = makeoutputfolder(Model,Runtype)

%% [outputdir,outputdir0,runID] = makeoutputfolder(Model,Runtype)
% Generate outputfolders for results. Outputdir0 is directory corresponding
% to current model settings, while outputdir also includes runID. 

% Mai Winstrup
% 2014-06-17 21:24: First independent version
% 2014-08-14 13:32: modelnumber -> runnumber
% 2014-08-19 16:18: Changes to runnumber, introduced runnumber=0, and 
%                   remove content of folder.
% 2014-08-21 20:01: Also works for SyntheticData, new name
% 2014-08-23 12:48: Removed changes to Model for synthetic data
% 2014-09-29 19:00: runnumber -> runID, outputdir0 for syn data changed.
% 2014-10-07 14:57: Nolonger saving runID=0

%% Find/generate output folder:
% Running in development mode?
if strcmp(Runtype.develop,'yes'); outputdir0 = './Output/develop';
else outputdir0 = './Output';
end

if (strcmp(Model.icecore, 'SyntheticData'))
        outputdir0 = [outputdir0 '/SyntheticData/' Model.info];
else       
        all_species = Model.species{1};
        for j = 2:Model.nSpecies
            all_species = [all_species '_' Model.species{j}];
        end
        outputdir0 = [outputdir0 '/' Model.icecore '/Timescale/' ...
            num2str(Model.dstart) '-' num2str(Model.dend) 'm/' all_species];
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