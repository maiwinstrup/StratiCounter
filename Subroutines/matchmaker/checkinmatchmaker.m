function checkinmatchmaker(outputdir,Model,additionaldata,Runtype)

%% checkinmatchmaker(outputdir,Model,additionaldata,Runtype)
% Extracts timescale results; opens them in Matchmaker for manual check,
% where they are compared to the initial set of manual counts (if existing). 
% The function can also be used to compare the results of several runs of 
% the algorithm. Outputdir is then a cell array with multiple cells.  
% Manual counts (if such exists) are provided with the background of all 
% data series available. If wanted to include additional data series (e.g.
% sulfate) to the output counts, this can be done by passing the names of 
% these as cell array "additionaldata". 

% Copyright (C) 2015  Mai Winstrup

%% Convert outputdir to cell array if required: 
if ischar(outputdir)
    outputdir1{1}=outputdir;
    outputdir = outputdir1;
end

%% Folder for matchfiles:
if strcmp(Runtype.develop,'yes')
    matchfilefolder = ['./Output/develop/' Model.icecore '/Matchfiles'];
else
    matchfilefolder = ['./Output/' Model.icecore '/Matchfiles'];
end
if ~exist(matchfilefolder,'dir'); mkdir(matchfilefolder); end

%% 1a: Manual layer counts
% Convert layer counts to format used in matchmaker:
manualcounts = loadlayercounts(Model,[Model.dstart Model.dend]);

% Must have more than 5 manual counts for these to be displayed:
if size(manualcounts,1)>5
    mp(:,1) = manualcounts(:,1);
    mp(:,2) = 1; % Certain years
    mask = manualcounts(:,3)==1;
    mp(mask,2) = 2; % Uncertain years
    
    % Save layer counts:
    filename = [matchfilefolder '/layers_manual.mat'];
    save(filename,'mp')

    %% 1b: Matchmaker file including all available data files:
    load(Model.path2data)
    allspecies = data.name; % All species are to be included.

    % Number of species for core
    nSp(1) = length(allspecies);

    % Save data file in matchmaker format:
    matchmakerdata(allspecies,Model,'_all');
    clear data

    %% 1c: Start matchmaker's files_main:
    % Create files_main file
    fid = fopen('./Subroutines/matchmaker/files_main.m','w'); % open new file, discard any content

    % Manual counts with all available data records:
    iPanel0 = 1;
    textinfile = ['files.core{' num2str(iPanel0) '}=''' Model.nameManualCounts(1:end-4) '''; \r\n'...
        'files.datafile{' num2str(iPanel0) '}=''data_all.mat''; \r\n'...
        'files.matchfile{' num2str(iPanel0) '}=''' matchfilefolder '/layers_manual.mat''; \r\n \r\n'];
    fprintf(fid,textinfile);
    fclose(fid);

else
    % Create files_main file without content:
    fid = fopen('./Subroutines/matchmaker/files_main.m','w'); % open new file, discard any content
    fclose(fid);
    iPanel0 = 0;
end

%% 2: Output of algorithm:
% Number of output layerings: 
nOutput = length(outputdir);
% Number of species displayed: 
nSp(iPanel0+1:iPanel0+nOutput) = nan;

for iOutput = 1:nOutput
    iPanel = iOutput+iPanel0;
    clear mp data species depth depth_no colours

    %% Load timescale and model:
    load([outputdir{iOutput} '/Model.mat'])
    timescale = importdata([outputdir{iOutput} '/' Model.icecore ...
        '_timescale_' num2str(Model.dstart) '-' num2str(Model.dend) 'm.txt']);
    timescale1yr = timescale.data;
   
    % Convert layer counts to format used in matchmaker:
    mp(:,1) = timescale1yr(:,1);
    mp(:,2) = 1; % Certain layers, large grey bars

    %% Tiepoints and marker horizons:
    % Tiepoints are added too:  
    tp1 = Model.tiepoints; % actual tiepoints, used in derivation
    if ~isempty(tp1)
        mp = [mp; tp1(:,1), ones(length(tp1),1)*5]; % Small blue bars
    end

    % Marker horizons:
    if ~isempty(Model.dMarker)
        for i = 1:length(Model.dMarker)
            tp = Model.dMarker{i};
            mp = [mp; tp(:), ones(length(tp),1)*5]; % Small blue bars
        end
    end
    
    % Sort:
    mp = sortrows(mp);

    % Save layer counts and tiepoints:
    if iOutput == 1
        filename = [matchfilefolder '/layers_auto.mat']; 
    else
        filename = [matchfilefolder '/layers_auto' num2str(iOutput) '.mat'];
    end
    save(filename,'mp')
    
    %% Make matchmaker datafile (unprocessed data):
    % Add additional data series, if so desired:
    Model.species=[Model.species additionaldata];
    Model.nSpecies = length(Model.species);
    if iOutput == 1
        matchmakerdata(Model.species,Model,'_auto')
    else
        matchmakerdata(Model.species,Model,['_auto' num2str(iOutput)])
    end
    
    %% Create files_main file:
    fid = fopen('./Subroutines/matchmaker/files_main.m','a'); % append to file

    % Convert to text:
    if iOutput == 1
        textinfile = ['files.core{' num2str(iPanel) '}=''' Model.icecore '''; \r\n'...
            'files.datafile{' num2str(iPanel) '}=''data_auto.mat''; \r\n'...
            'files.matchfile{' num2str(iPanel) '}=''' matchfilefolder '/layers_auto.mat''; \r\n \r\n'];
    else
        textinfile = ['files.core{' num2str(iPanel) '}=''' Model.icecore '''; \r\n'...
            'files.datafile{' num2str(iPanel) '}=''data_auto' num2str(iOutput) '.mat''; \r\n'...
            'files.matchfile{' num2str(iPanel) '}=''' matchfilefolder '/layers_auto' num2str(iOutput) '.mat''; \r\n \r\n'];
    end
    fprintf(fid,textinfile);
    fclose(fid);
    
    % Nspecies for core
    nSp(iPanel) = Model.nSpecies+length(additionaldata);
end

%% Open matchmaker:
% Total number of panels displayed: 
nPanel = nOutput+iPanel0;

% Maximum number of data records displayed per panel:
nSpMax = round(10/nPanel);
nSp = min(nSp,nSpMax);

%% Generate matchmaker_sett.mat:
% Set x-limits:
sett.xlim(1:2,1)=Model.dstart;
if ~isempty(manualcounts)
    meanlambda = mean(diff(manualcounts(:,1)));
else
    % Estimated based on last set of autocounts:
    meanlambda = mean(diff(mp(:,1))); 
end

if Model.dend < Model.dstart+30*meanlambda
    sett.xlim(1:2,2) = Model.dend;
else
    sett.xlim(1:2,2)=Model.dstart+round(20*meanlambda);
end
sett.specs{1} = [1];
sett.specs{2} = [1];
save('./Subroutines/matchmaker/matchmaker_sett.mat','sett','-v6')

% Open matchmaker:
matchmaker('files_main',1:nPanel,nSp)