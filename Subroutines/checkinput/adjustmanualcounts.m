function [manualcounts, Model] = adjustmanualcounts(manualcounts,Model)

%% [manualcounts, Model] = adjustmanualcounts(manualcounts,Model)
% Check and adjust the format of the manual layer counts. Ages are 
% converted to ageUnitOut. 
% Copyright (C) 2015  Mai Winstrup

%% No changes if no layer counts are given:
if isempty(manualcounts);
    return; 
end

%% Ensure that the increase/decreasing age numbers correspond to the
% selected value of ageUnitManual:
if manualcounts(2,2)-manualcounts(1,2)>0 
    % An increase in age values with increasing depth is not consistent 
    % with manual counts age unit "AD"
    if strcmp(Model.ageUnitManual,'AD')
        disp(['OBS: Age values are increasing with depth, which is '...
            'inconsistent with the value of Model.ageUnitManual ("AD"):'])
        disp('     Model.ageUnitManual is changed to "layers"');
        Model.ageUnitManual = 'layers';
    end
end

%% Uncertainties of the individual layers: 
% Must either be ones or zeros:
mask = manualcounts(:,3)~=0 & manualcounts(:,3)~=1;
if sum(mask)>0
    disp(['Check manual layer uncertainties. These must either be given a ' ...
        'value of 0 (certain) or 1 (uncertain)']);
    return
end

%% Accumulated uncertainty within current depth interval: 
% Must always start with an uncertainty of 0:
% First certain layer: 
first_certain_layer = find(manualcounts(:,3)==0,1,'first');
manualcounts(:,4)=manualcounts(:,4)-manualcounts(first_certain_layer,4);

%% Converting manual layer counts to "ageUnitOut":
% First convert to the b2k age convention. This is done unless only the 
% layer number relative to the chosen start depth is given, in which case 
% this is, of course, not possible.
switch Model.ageUnitManual
    case 'BP'
        manualcounts(:,2) = manualcounts(:,2)+50;
    case 'AD'
        manualcounts(:,2) = [2000-manualcounts(2:end,2)-1; 2000-manualcounts(end,2)];        
    case 'layers'
        % Counting the layers from the chosen start depth:
        manualcounts(:,2) = manualcounts(:,2)-manualcounts(first_certain_layer,2)+1;
    case 'b2k'
    otherwise
        disp('ageUnitManual not recognized')
end

% If only relative ages are given, age convention of the output must also
% be relative:
% If tiepoints exist, these are used to provide ages:
if ~isempty(Model.tiepoints) 
    if strcmp(Model.ageUnitTiepoints,'layers')
        Model.ageUnitOut = 'layers';
        disp('Age units of resulting timescale is changed to "layers"')
    end
end
if strcmp(Model.ageUnitManual,'layers')
    Model.ageUnitOut = 'layers';
    disp('Age units of resulting timescale is changed to "layers"')
end
        
%% Then convert to ageUnitOut:
switch Model.ageUnitOut
    case 'AD'
        manualcounts(:,2) = [2000-manualcounts(1,2); 2000-manualcounts(1:end-1,2)-1];
    case 'BP'
        manualcounts(:,2) = manualcounts(:,2)-50;
    case 'layers'
        manualcounts(:,2) = manualcounts(:,2)-manualcounts(first_certain_layer,2)+1;
    case 'b2k'
    otherwise
        disp('ageUnitOut not recognized')
end

%% If we have uncertain years, manually counted ages may start at "x.5 yr". 
% We do not wish our new timescale to start with a half year (all years 
% would then be half years).
% In this case, all ages are (subtracted) by 0.5. 
% Examples: 5836.5b2k -> 5836 b2k; 
%           -3836.5AD -> -3837AD; 
%           277.5AD -> 277 AD
layerfraction = mod(manualcounts(first_certain_layer,2),1);
if layerfraction ~= 0
    manualcounts(:,2) = manualcounts(:,2)-layerfraction;
    disp(['Manual ages is subtracted by ' num2str(layerfraction) ' year '...
        'in order to provide integer years'])
end