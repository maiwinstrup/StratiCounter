function checkinmatchmaker(outputdir,Model,Runtype)

%% checkinmatchmaker(outputdir,Model)
% Extract timescale results; open them in Matchmaker for manual check,
% and compare to initial manual counts.
% Copyright (C) 2015  Mai Winstrup

%% 1a: Manual layer counts
% Convert layer counts to format used in matchmaker:
manualcounts = loadlayercounts(Model,[Model.dstart Model.dend]);

mp(:,1) = manualcounts(:,1);
mp(:,2) = 1; % Certain years
mask = manualcounts(:,3)==1;
mp(mask,2) = 2; % Uncertain years

% Save layer counts:
if strcmp(Runtype.develop,'yes')
    matchfilefolder = ['./Output/develop/' Model.icecore '/Matchfiles'];
else
    matchfilefolder = ['./Output/' Model.icecore '/Matchfiles'];
end
if ~exist(matchfilefolder,'dir'); mkdir(matchfilefolder); end
filename = [matchfilefolder '/layers_manual.mat'];
save(filename,'mp')

%% 1b: Matchmaker file including all available data files:
load(Model.path2data)
allspecies = data.name; % All species are to be included.

% Save data file in matchmaker format:
matchmakerdata(allspecies,Model,'_manual');
clear data

%% 1c: Start matchmaker's files_main:
% Create files_main file
fid = fopen('./Subroutines/matchmaker/files_main.m','w'); % open new file, discard any content

% Convert to text:
iCore = 1;
textinfile = ['files.core{' num2str(iCore) '}=''' Model.nameManualCounts(1:end-4) '''; \r\n'...
    'files.datafile{' num2str(iCore) '}=''data_manual.mat''; \r\n'...
    'files.matchfile{' num2str(iCore) '}=''' matchfilefolder '/layers_manual.mat'';'];
fprintf(fid,textinfile);
fclose(fid);
    
% Nspecies for core
nSp(iCore) = length(allspecies);

%% 2: Output of algorithm:
dir{2} = outputdir; % Current 
% We may add others too:
% dir{3} = ... 

% Sometimes we may want to include additional data series, e.g. sulfate:
additionaldata{2} = [];

nCore = length(dir);
for iCore = 2:nCore
    clear mp data species depth depth_no colours

    %% Load timescale and model:
    load([dir{iCore} '/Model.mat'])
    timescale = importdata([dir{iCore} '/' Model.icecore '_timescale.txt']);
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
    if iCore == 2
        filename = [matchfilefolder '/layers_auto.mat']; ...
            % We may compare two results from the same ice core
    else
        filename = [matchfilefolder '/layers_auto' num2str(iCore-1) '.mat'];...
    end
    save(filename,'mp')
    
    %% Make matchmaker datafile (unprocessed data):
    % Add additional data series, if so desired:
    Model.species=[Model.species additionaldata{iCore}];
    Model.nSpecies = length(Model.species);
    if iCore == 2
        matchmakerdata(Model.species,Model,'_auto')
    else
        matchmakerdata(Model.species,Model,['_auto' num2str(iCore-1)])
    end
    
    %% Create files_main file:
    fid = fopen('./Subroutines/matchmaker/files_main.m','a'); % append to file

    % Convert to text:
    if iCore == 2
        textinfile = ['\r\nfiles.core{' num2str(iCore) '}=''' Model.icecore '''; \r\n'...
            'files.datafile{' num2str(iCore) '}=''data_auto.mat''; \r\n'...
            'files.matchfile{' num2str(iCore) '}=''' matchfilefolder '/layers_auto.mat''; \r\n'];
    else
        textinfile = ['\r\nfiles.core{' num2str(iCore) '}=''' Model.icecore '''; \r\n'...
            'files.datafile{' num2str(iCore) '}=''data_auto' num2str(iCore-1) '.mat''; \r\n'...
            'files.matchfile{' num2str(iCore) '}=''' matchfilefolder '/layers_auto' num2str(iCore-1) '.mat''; \r\n'];
    end
    fprintf(fid,textinfile);
    fclose(fid);
    
    % Nspecies for core
    nSp(iCore) = Model.nSpecies;
end

%% Open matchmaker:
% Display a maximum of 5 subfigures per panel:
nSp = min(nSp,5);

% And a total number of maximum 10:
while sum(nSp) > 10
    [~,index]=max(nSp);
    nSp(index) = nSp(index)-1;
end

%% Generate matchmaker_sett.mat:
% Set x-limits:
sett.xlim(1:2,1)=Model.dstart;
meanlambda = mean(diff(manualcounts(:,1)));
if Model.dend < Model.dstart+30*meanlambda
    sett.xlim(1:2,2) = Model.dend;
else
    sett.xlim(1:2,2)=Model.dstart+round(20*meanlambda);
end
sett.specs{1} = [1];
sett.specs{2} = [1];
save('./Subroutines/matchmaker/matchmaker_sett.mat','sett','-v6')

% Open matchmaker:
matchmaker('files_main',1:nCore,nSp)