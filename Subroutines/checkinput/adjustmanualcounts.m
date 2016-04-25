function [manualcounts, Model] = adjustmanualcounts(manualcounts,Model)

%% [manualcounts, Model] = adjustmanualcounts(manualcounts,Model)
% Check and adjust the format of the manual layer counts. Ages are 
% converted to Model.Out.ageUnit. 
% Copyright (C) 2015  Mai Winstrup

%% No changes if no layer counts are given:
if isempty(manualcounts);
    return; 
end

%% Ensure that increasing/decreasing age numbers are consistent with the 
% value of ageUnitManual:
if size(manualcounts,1)>1 
    if manualcounts(2,2)-manualcounts(1,2)>0 
        % An increase in age values with increasing depth is not consistent 
        % with manual counts age unit "AD"
        if strcmp(Model.ageUnitManual,'AD')
            warning(['Manual timescale is inconsistent with the value of '...
                'Model.ageUnitManual ("AD"). '...
                'Model.ageUnitManual is changed to "layers"']);
            Model.ageUnitManual = 'layers';
        end
    elseif manualcounts(2,2)-manualcounts(1,2)<0
        if ~strcmp(Model.ageUnitManual,'AD')
            warning(['Manual timescale is inconsistent with the value of '... 
                'Model.ageUnitManual. '...
                'Model.ageUnitManual is changed to "AD"']);
            Model.ageUnitManual = 'AD';
        end
    end
end

%% Uncertainties of the individual layers: 
% Must either be ones or zeros:
mask = manualcounts(:,3)~=0 & manualcounts(:,3)~=1;
if sum(mask)>0
    warning(['Check manual layer uncertainties. These must either be given '...
        'a value of 0 (certain) or 1 (uncertain). For now, they have been '...
        'assumed certain.']);
    manualcounts(mask,:) = 0;
end

%% Accumulated uncertainty within current depth interval: 
% Must always start with an uncertainty of 0:
% First certain layer: 
first_certain_layer = find(manualcounts(:,3)==0,1,'first');
manualcounts(:,4)=manualcounts(:,4)-manualcounts(first_certain_layer,4);

%% Converting units of manual layer counts to "ageUnitOut":

%% Using tiepoints to convert from relative manual ages:
% When needed (and possible), tiepoints are used to convert relative manual 
% ages to absolute manual ages:
if ~isempty(Model.tiepoints) 
    if strcmp(Model.ageUnitManual,'layers')
        if strcmp(Model.ageUnitTiepoints,'layers') 
            if ~strcmp(Model.Out.ageUnit,'layers')
                Model.Out.ageUnit = 'layers';
                disp('Age units of resulting timescale is changed to "layers"')
            end
            
        else
            % Convert the relative ages of the manual layer counts to 
            % absolute ages by fitting to the first tiepoint.             
            % We extend manual layer counts by adding an extra layer at the 
            % beginning and end of manual counts: 
            lambda_before = 1.5*(manualcounts(2,1)-manualcounts(1,1));
            lambda_after = 1.5*(manualcounts(end,1)-manualcounts(end-1,1));
            switch Model.ageUnitManual
                case 'AD'
                    manualcounts_add = [manualcounts(1,1)-lambda_before, ...
                        manualcounts(1,2)+1, 0, 0; manualcounts; ...
                        manualcounts(end,1)+lambda_after, manualcounts(end,2)-1, 0, 0];
                otherwise 
                    manualcounts_add = [manualcounts(1,1)-lambda_before, ...
                        manualcounts(1,2)-1, 0, 0; manualcounts; ...
                        manualcounts(end,1)+lambda_after, manualcounts(end,2)+1, 0, 0];
            end
            % Layer number for first tiepoint:
            tplayernumber = floor(interp1(manualcounts_add(:,1),...
                manualcounts_add(:,2),Model.tiepoints(1,1))); 
            clear manualcounts_add
            % This value will be nan if manual layers do not cover the
            % first tiepoint. 
            
            % Convert relative ages to absolute values (if possible):
            if isfinite(tplayernumber)
                switch Model.ageUnitTiepoints
                    case 'AD'
                        startage = Model.tiepoints(1,2)+(tplayernumber-manualcounts(1,2)+1);
                        manualcounts(:,2) = startage-(manualcounts(:,2)-manualcounts(1,2)); 
                    otherwise 
                        startage = Model.tiepoints(1,2)-(tplayernumber-manualcounts(1,2))-1;
                        manualcounts(:,2) = startage+(manualcounts(:,2)-manualcounts(1,2)); 
                end
                % Record this change in age unit for the manual timescale: 
                Model.ageUnitManual = Model.ageUnitTiepoints;
            else
                if ~strcmp(Model.Out.ageUnit,'layers')
                    Model.Out.ageUnit = 'layers';
                    disp('Age units of resulting timescale is changed to "layers"')
                end
            end
        end
    end
end

%% Convert manual ages to the b2k age convention. 
% This is done unless only the layer number relative to the chosen start 
% depth is given, in which case this is, of course, not possible.
switch Model.ageUnitManual
    case 'BP'
        manualcounts(:,2) = manualcounts(:,2)+50;
    case 'AD'
        manualcounts(:,2) = [2000-manualcounts(2:end,2)-1; 2000-manualcounts(end,2)];        
    case 'layers'
        % Counting the layers from the chosen start depth:
        if isempty(Model.tiepoints)
            manualcounts(:,2) = manualcounts(:,2)-manualcounts(first_certain_layer,2)+1;
        else
            manualcounts(:,2) = manualcounts(:,2)-manualcounts(1,2)+1;
        end
    case 'b2k'
    otherwise
        disp('ageUnitManual not recognized')
end

%% If possible, relative tiepoint ages are converted to absolute 
% ages:
if ~isempty(Model.tiepoints) 
    if strcmp(Model.ageUnitTiepoints,'layers')
        if ~strcmp(Model.ageUnitManual,'layers')
            % Convert the relative ages of the tiepoints to absolute ages 
            % by fitting to the manual timescale (now in b2k). 
            
            % Extend manual layer counts (in b2k) by adding an extra layer 
            % at the beginning and end of manual counts: 
            lambda_before = 1.5*(manualcounts(2,1)-manualcounts(1,1));
            lambda_after = 1.5*(manualcounts(end,1)-manualcounts(end-1,1));
            manualcounts_add = [manualcounts(1,1)-lambda_before, ...
                manualcounts(1,2)+1, 0, 0; manualcounts; ...
                manualcounts(end,1)+lambda_after, manualcounts(end,2)-1, 0, 0];
            
            % Ages (in b2k, no rounding) at tiepoints:
            tpages = interp1(manualcounts_add(:,1),...
                manualcounts_add(:,2),Model.tiepoints(:,1)); 
            if isfinite(tpages)
                Model.tiepoints(:,2) = tpages;
            elseif ~strcmp(Model.Out.ageUnit,'layers')
                Model.Out.ageUnit = 'layers';
                disp('Age units of resulting timescale is changed to "layers"')
            end
        end
    else
        % Convert tiepoint ages to the b2k age convention:
        switch Model.ageUnitTiepoints
            case 'BP'
                Model.tiepoints(:,2) = Model.tiepoints(:,2)+50;
            case 'AD'
                Model.tiepoints(:,2) = 2000-Model.tiepoints(:,2);
        end
    end
end

%% Then convert to ageUnitOut:
% Manual layer counts:
switch Model.Out.ageUnit
    case 'AD'
        manualcounts(:,2) = [2000-manualcounts(1,2); 2000-manualcounts(1:end-1,2)-1];
    case 'BP'
        manualcounts(:,2) = manualcounts(:,2)-50;
    case 'layers'
        manualcounts(:,2) = manualcounts(:,2)-manualcounts(first_certain_layer,2)+1;
    case 'b2k'
    otherwise
        disp('Model.Out.ageUnit not recognized')
end

% Tiepoints:
if ~isempty(Model.tiepoints)
    switch Model.Out.ageUnit
        case 'AD'
            Model.tiepoints(:,2) = 2000-Model.tiepoints(:,2);
        case 'BP'
            Model.tiepoints(:,2) = Model.tiepoints(:,2)-50;
        case 'layers'
            Model.tiepoints(:,2) = Model.tiepoints(:,2)-Model.tiepoints(1,2)+1;
    end
    % Round tiepoint ages:
    Model.tiepoints(:,2) = floor(Model.tiepoints(:,2));
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