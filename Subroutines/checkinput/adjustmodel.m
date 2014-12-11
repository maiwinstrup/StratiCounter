function Model = adjustmodel(Model)

%% Model = adjustmodel(Model)
% This function checks that the Model structure array has the correct 
% format, and (if possible) makes the required corrections to the array. 
% Mai Winstrup, 2014
% 2014-06-16 12:09

%% Adjust depth interval for layer counting:
% The depths below do not include displacement due to possible non-zero
% value of dx_center:
Model.dstart = ceil((Model.dstart-Model.dx_center)/Model.dx)*Model.dx; 
Model.dend = ceil((Model.dend-Model.dx_center)/Model.dx)*Model.dx; 

%% Same for manual templates:
Model.manualtemplates = ceil((Model.manualtemplates-Model.dx_center)/Model.dx)*Model.dx; 


%% Ensure that Model.dstart and Model.dend corresponds to a data point: nej
%Model.manualtemplates = checkdepth(Model.manualtemplates,Model.dx,Model.dx_center);
Model.initialpar = checkdepth(Model.initialpar,Model.dx,Model.dx_center);

%% Model initialpar : bør være indenfor data interval. Der er ikke nogen 
% grund til andet, da det alligevel kun er estimat, som ændres senere. 

%%
Model.deriv = sort(Model.deriv);
% Model.deriv bør altid indeholder data serien selv. 
% (forudsættes f.s i finalizepreprocessing) og det giver mening!!

%% Sort ordering of data species:
[Model.species, index] = sort(Model.species);
% And corresponding preprocessing steps etc.:
Model.preprocess = Model.preprocess(index);
Model.wSpecies = Model.wSpecies(index);
clear index

%% Partition the preprocessing steps into initial and recurring processing 
% steps (i.e. those with floating distances). 
for j = 1:Model.nSpecies
    if size(Model.preprocess{j},2)>2
        disp('check preprocessing format!')
        return
    end
end

index = zeros(1,Model.nSpecies); % Last preprocessing step with fixed distance
for j = 1:Model.nSpecies
    stringpos = strfind(Model.preprocess{j}(:,1),'float');
    for istep = 1:size(Model.preprocess{j},1)
        if isempty(stringpos{istep}); index(j) = istep; 
        else break
        end
    end
    Model.preprocess_init{j} = Model.preprocess{j}(1:index(j),:);
    Model.preprocess_rep{j} = Model.preprocess{j}(index(j)+1:end,:);
   % clear Model.preprocess
   
    % Replace with "none" if array is empty:
    if isempty(Model.preprocess_init{j})
        Model.preprocess_init{j}= {'none',[]};
    end
    if isempty(Model.preprocess_rep{j})
        Model.preprocess_rep{j}= {'none',[]};
    end
end

%% Does one data species occur more than once?
nUnique = length(unique(Model.species));
if nUnique~=Model.nSpecies
    disp('One data species occur more than once, please correct')
    return
end

%% Model.wSpecies must be the correct format:
Model.wSpecies = Model.wSpecies(:);
% And of correct length:
if length(Model.wSpecies) ~= Model.nSpecies
    disp('Model.wSpecies does not have the required format, please correct')
    return
end

%% Depth interval must be positive:
if Model.dend <= Model.dstart
    disp('Check depth interval')
    return
end
% Model.dstart skal være et data punkt?

%% If tiepoints: 
if ~isempty(Model.tiepoints)
    % Tiepoints are sorted relative to depth:
    Model.tiepoints = sortrows(Model.tiepoints,1);

    % Remove tiepoints from outside data interval:
    mask = Model.tiepoints(:,1)>=Model.dstart & Model.tiepoints(:,1)<=Model.dend;
    Model.tiepoints = Model.tiepoints(mask,:);
    
    % Only the section between the uppermost and lowermost tiepoints 
    % (within given depth interval) is considered. Values of dstart and 
    % dend are changed to reflect this.
    Model.dstart = floor(Model.tiepoints(1,1));
    Model.dend = ceil(Model.tiepoints(end,1));
    
    % Check the age unit of tiepoints:
    while ~ismember(Model.ageUnitTiepoints,{'AD','BP','b2k','layers'})
        promt = ['Which timescale terminology was used for tiepoints?' ... 
            '\n(Options: AD, BP, b2k, layers): '];
        Model.ageUnitTiepoints = input(promt,'s');
    end

    % Convert tiepoints to integer values:
%    switch Model.ageUnitTiepoints
        Model.tiepoints(:,2)=floor(Model.tiepoints(:,2));        
%        case 'AD'; Model.tiepoints(:,2)=floor(Model.tiepoints(:,2));        
%        otherwise; Model.tiepoints(:,2)=ceil(Model.tiepoints(:,2));
%    end
end

%% If using 'FFT' as Model.type: 
% Value of model order should be increased. 
% In this case, the original value denotes the number of phase-shifted 
% sinusoides, thus giving rise to a total number of parameters equal to: 
if strcmp(Model.type,'FFT')
    Model.order = 1+2*Model.order;
end

%% Check depth interval for layer characteristics:
if Model.manualtemplates(1) >= Model.manualtemplates(2)
    disp('Check depth interval for layer templates')
end
% And for initial layer parameters:
if Model.initialpar(1) >= Model.initialpar(2)
    disp('Check depth interval for initial layer parameter estimation')
end

%% Sections for mean layer thickness calculations:
% Remove sections from outside depth interval:
for i = 1:length(Model.dMarker)
    mask = Model.dMarker{i}>=Model.dstart & Model.dMarker{i}<=Model.dend;
    Model.dMarker{i} = Model.dMarker{i}(mask);
end

%% Check that one of four options are chosen as unit for timescale output: 
while ~ismember(Model.ageUnitOut,{'AD','BP','b2k','layers'})
    promt = 'Which timescale terminology to be used for output? \n(Options: AD, BP, b2k, layers): ';
    Model.ageUnitOutput = input(promt,'s');
end

%% Truncating tiepoints etc. to the desired data resolution:
% Tiepoints:
if ~isempty(Model.dx)
    depth_new = Model.dstart+Model.dx_center:Model.dx:Model.dend;
    if ~isempty(Model.tiepoints)
        % Interpolate to data resolution:
        index = interp1(depth_new,1:length(depth_new),Model.tiepoints(:,1),'nearest','extrap');
        Model.tiepoints(:,1) = depth_new(index);
    end

    % Same for sections for lambda calculations:
    for i = 1:length(Model.dMarker)
        mask = Model.dMarker{i}>=Model.dstart & Model.dMarker{i}<=Model.dend;
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

%% Check også, at vi ikke har f.x. process = {minusmin,2] eller [log, 2];
