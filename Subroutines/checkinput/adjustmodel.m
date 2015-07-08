function Model = adjustmodel(Model)

%% Model = adjustmodel(Model)
% This function checks that the Model structure array has the correct 
% format, and (if possible) makes the required corrections to the array. 
% Copyright (C) 2015  Mai Winstrup

%% Check value of dx_center:
if Model.dx_center<0 || Model.dx_center>=1
    disp('Value of Model.dx_center must be a positive number below 1. Please correct!')
    return
end
if Model.dx_center~=0 && Model.dx_center~=0.5;
    disp(['Check value of Model.dx_center (current value: ' num2str(Model.dx_center) ')'])
end

%% Adjust depth intervals:
% The depths below do not include displacement due to possible non-zero
% value of dx_center.
% Interval for layer counting:
Model.dstart = ceil((Model.dstart-Model.dx_center*Model.dx)/Model.dx)*Model.dx; 
Model.dend = ceil((Model.dend-Model.dx_center*Model.dx)/Model.dx)*Model.dx; 

% Interval for manual templates:
Model.manualtemplates = ...
    ceil((Model.manualtemplates-Model.dx_center*Model.dx)/Model.dx)*Model.dx; 

% Interval for initial parameters:
Model.initialpar = ...
    ceil((Model.initialpar-Model.dx_center*Model.dx)/Model.dx)*Model.dx; 

% Depth intervals must be positive:
if Model.dend<=Model.dstart
    disp('Check depth interval for layer counting')
    return
end
if Model.manualtemplates(1)>=Model.manualtemplates(2)
    disp('Check depth interval for layer templates')
end
if Model.initialpar(1)>=Model.initialpar(2)
    disp('Check depth interval for initial layer parameter estimation')
end

%% Format preprocessing steps:
for j = 1:Model.nSpecies
    for k = 1:2
        if size(Model.preprocsteps{j,k},2)==1
            Model.preprocsteps{j,k} = {Model.preprocsteps{j,k}{1,1},[]};
        end
    end
end

%% Sort ordering of data species:
[Model.species, index] = sort(Model.species);

% Similar for corresponding preprocessing steps and weights:
Model.preprocsteps = Model.preprocsteps(index,:);
Model.wSpecies = Model.wSpecies(index);
clear index

%% Does one data species occur more than once?
nUnique = length(unique(Model.species));
if nUnique~=Model.nSpecies
    disp('One data species occur more than once, please correct')
    return
end

%% Model.wSpecies must be the correct format:
Model.wSpecies = Model.wSpecies(:);
% And of correct length:
if length(Model.wSpecies)~=Model.nSpecies
    disp('Model.wSpecies does not have the required format, please correct')
    return
end

%% Tiepoints: 
if ~isempty(Model.tiepoints)
    % Tiepoints are sorted relative to depth:
    Model.tiepoints = sortrows(Model.tiepoints,1);

    % Remove tiepoints from outside data interval:
    mask = Model.tiepoints(:,1)>=Model.dstart+Model.dx_center*Model.dx &...
        Model.tiepoints(:,1)<=Model.dend+Model.dx_center*Model.dx;
    Model.tiepoints = Model.tiepoints(mask,:);
    
    % Only the section between the uppermost and lowermost tiepoints 
    % (within given depth interval) is considered. Values of dstart and 
    % dend are changed to reflect this.
    Model.dstart = floor(Model.tiepoints(1,1));
    Model.dend = ceil(Model.tiepoints(end,1));
    
    % Check age unit of tiepoints:
    if ~isfield(Model,'ageUnitTiepoints')
        promt = ['Which timescale terminology was used for tiepoints?' ... 
            '\n(Options: AD, BP, b2k, layers): '];
        Model.ageUnitTiepoints = input(promt,'s');
    end        
    while ~ismember(Model.ageUnitTiepoints,{'AD','BP','b2k','layers'})
        promt = ['Which timescale terminology was used for tiepoints?' ... 
            '\n(Options: AD, BP, b2k, layers): '];
        Model.ageUnitTiepoints = input(promt,'s');
    end

    % Convert tiepoints to integer values:
    Model.tiepoints(:,2)=floor(Model.tiepoints(:,2));        
end

%% If using 'FFT' as Model.type: 
% Value of model order should be increased. 
% In this case, the original value denotes the number of phase-shifted 
% sinusoides, thus giving rise to a total number of parameters equal to: 
if strcmp(Model.type,'FFT')
    Model.order = 1+2*Model.order;
end

%% Sections for mean layer thickness calculations:
if ~iscell(Model.dMarker)&&~isempty(Model.dMarker)
    Model.dMarker = {Model.dMarker};
end

% Remove sections from outside depth interval:
for i = 1:length(Model.dMarker)
    mask = Model.dMarker{i}>=Model.dstart & Model.dMarker{i}<=Model.dend;
    Model.dMarker{i} = Model.dMarker{i}(mask);
end

%% Check for file extension on filename with manual layer counts:
if ~strcmp(Model.nameManualCounts(end-3:end),'.txt')
    % If not: add to filename
    Model.nameManualCounts = [Model.nameManualCounts '.txt'];
end

%% Check that one of four options are chosen as unit for timescale output: 
while ~ismember(Model.ageUnitOut,{'AD','BP','b2k','layers'})
    promt = 'Which timescale terminology to be used for output? \n(Options: AD, BP, b2k, layers): ';
    Model.ageUnitOutput = input(promt,'s');
end

%% Truncating tiepoints etc. to the desired data resolution:
if ~isempty(Model.dx)
    % New depth scale:
    depth_new = Model.dstart+Model.dx_center*Model.dx:Model.dx:...
        Model.dend+Model.dx_center*Model.dx;
    
    % Tiepoints:
    if ~isempty(Model.tiepoints)
        % Interpolate to data resolution:
        index = interp1(depth_new,1:length(depth_new),Model.tiepoints(:,1),'nearest','extrap');
        Model.tiepoints(:,1) = depth_new(index);
    end

    % Sections for lambda calculations:
    for i = 1:length(Model.dMarker)
        mask = Model.dMarker{i}>=Model.dstart+Model.dx_center*Model.dx & ...
            Model.dMarker{i}<=Model.dend+Model.dx_center*Model.dx;
        index = interp1(depth_new,1:length(depth_new),Model.dMarker{i}(mask),'nearest');
        Model.dMarker{i} = depth_new(index);
    end
end

if ~isempty(Model.dMarker)
    if isempty(Model.dMarker{1})
        Model.dMarker = [];
    else
        for i = 1:length(Model.dMarker)
            Model.dMarker{i} = sort(Model.dMarker{i});
        end
    end
end

%% Confidence interval:
for i = 1:Model.confInterval
    quantile_perc = [(100-Model.confInterval)/2, 100-(100-Model.confInterval)/2];
    quantile_perc = sort(quantile_perc);
end
Model.prctile = quantile_perc/100;