function checkinmatchmaker(outputdir,Model)

%% checkinmatchmaker(outputdir,Model)
% Extract timescale results; open them in Matchmaker for manual check,
% and compare to initial manual counts.
% Mai Winstrup, 2014

%% 1a: Manual layer counts
% Convert layer counts to format used in matchmaker:
manualcounts = loadlayercounts(Model,[Model.dstart Model.dend]); % 2014-10-08 09:58

mp(:,1) = manualcounts(:,1);
mp(:,2) = 1; % Certain years
mask = manualcounts(:,2)==1;
mp(mask,2) = 2; % Uncertain years

% Save layer counts:
if ~exist('./matchfiles','dir'); mkdir('./matchfiles'); end
filename = ['./matchfiles/' Model.icecore 'layers_manual.mat'];
save(filename,'mp')

%% 1b: Matchmaker file including all available data files:
load(Model.pathData)
allspecies = data.name; % All species are to be included.

% Save data file in matchmaker format:
matchmakerdata(allspecies,Model,'manual');
clear data

%% 1c: Start matchmaker's files_main:
% Create files_main file
fid = fopen('files_main.m','w'); % open new file, discard any content

% Convert to text:
iCore = 1;
textinfile = ['files.core{' num2str(iCore) '}=''' 'manual''; \r\n'...
    'files.datafile{' num2str(iCore) '}=''' Model.icecore 'data_manual.mat''; \r\n'...
    'files.matchfile{' num2str(iCore) '}=''' Model.icecore 'layers_manual.mat'';'];
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
%additionaldata{2} = {'cond_Bern','Na_Bern','Ca_Bern','dust_Bern','NH4_Bern','NO3_Bern','H2O2_Bern','ECM','SO4','d18O','dD'}; 

nCore = length(dir);
for iCore = 2:nCore
    clear mp data species depth depth_no colours

    %% Load timescale and model:
    load([dir{iCore} '/Model.mat'])
    timescale = importdata([dir{iCore} '/' Model.icecore '_timescale1yr.txt']);
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
        filename = ['./matchfiles/' Model.icecore 'layers_auto_adj.mat']; % We may compare two results from the same ice core
    else
        filename = ['./matchfiles/' Model.icecore 'layers_auto' num2str(iCore-1) '_adj.mat']; % We may compare two results from the same ice core
    end
    save(filename,'mp')
    
    %% Make matchmaker datafile (unprocessed data):
    % Add additional data series, if so desired:
    Model.species=[Model.species additionaldata{iCore}];
    Model.nSpecies = length(Model.species);
    if iCore == 2
        matchmakerdata(Model.species,Model,'auto')
    else
        matchmakerdata(Model.species,Model,['auto' num2str(iCore-1)])
    end
    
    %% Create files_main file:
    fid = fopen('files_main.m','a'); % append to file

    % Convert to text:
    if iCore == 2
        textinfile = ['\r\nfiles.core{' num2str(iCore) '}=''' Model.icecore '''; \r\n'...
            'files.datafile{' num2str(iCore) '}=''' Model.icecore 'data_auto.mat''; \r\n'...
            'files.matchfile{' num2str(iCore) '}=''' Model.icecore 'layers_auto_adj.mat'';'];
    else
        textinfile = ['\r\nfiles.core{' num2str(iCore) '}=''' Model.icecore '''; \r\n'...
            'files.datafile{' num2str(iCore) '}=''' Model.icecore 'data_auto' num2str(iCore-1) '.mat''; \r\n'...
            'files.matchfile{' num2str(iCore) '}=''' Model.icecore 'layers_auto' num2str(iCore-1) '_adj.mat'';'];
    end
    fprintf(fid,textinfile);
    fclose(fid);
    
    % Nspecies for core
    nSp(iCore) = Model.nSpecies;
end

%% Open matchmaker:
% Set maximum 10 subfigures in total:
nSubfig = sum(nSp);
if nSubfig > 10
    nSp = round(min(nSp,10/nCore)); 
end

% Set x-limits:
sett.xlim(1:2,1)=Model.dstart;
meanlambda = mean(diff(manualcounts(:,1)));
sett.xlim(1:2,2)=Model.dstart+round(20*meanlambda);
sett.specs{1} = [1];
sett.specs{2} = [1];
save('matchmaker_sett.mat', 'sett', '-v6');

% Open matchmaker:
matchmaker('files_main',1:nCore,nSp)

end

function matchmakerdata(species,Model,datafilename)
%% matchmakerdata(species,Model,datafilename)
% Save unprocessed datafiles in matchmaker format. 

nSpecies = length(species);

data0 = loadrawdata(species{1},Model);
% Data:
data{1}=data0(:,2);
species{1}=species{1};
% Depth:
depth{1}=data0(:,1); 
depth_no(1)=1;
    
% Subsequent data files:    
for j = 2:nSpecies
    data0 = loadrawdata(species{j},Model);

    % Data: 
    data{j}=data0(:,2);
        
    % Depth scale:
    % On same depthscale as previous data files?
    for i = 1:length(depth)            
        if isequal(data0(:,1),depth{i})
            % Same depth scale:
            depth_no(j)=i;
            flag = 1; 
            break % Do not check further depth scales
        else
            flag = 0; 
        end
    end
    % Use old depthscale:
    if flag == 1
        continue % Continue to next data file
    end
     
    % New depth scale:
    nd = length(depth)+1; % Current "depth scale number"
    depth{nd}=data0(:,1);
    depth_no(j)=nd;
end
        
% Species are given random colours:
colours = rand(nSpecies,3);   
mask = colours>0.9;
colours(mask)=0.7*colours(mask); % avoiding too bleak colors
% Save data:
if ~exist('./data','dir'); mkdir('./data'); end
filename = ['./data/' Model.icecore 'data_' datafilename '.mat'];
save(filename,'data','depth','depth_no','species','colours')
end