function savemarkerastxt(markerConf,filename,Model)

%% savemarkerastxt(markerConf,Model,filename)
% Save most likely layer and confidence intervals of number of layers 
% between marker horizons as text file with header. 
% Copyright (C) 2015  Mai Winstrup

%% Meta data:
% Species and their associated weighting:
species_and_weight = [Model.species{1} ' (w=' num2str(Model.wSpecies(1)) ')'];
for j = 2:Model.nSpecies
    species_and_weight = [species_and_weight ', ' Model.species{j}  ...
        ' (w=' num2str(Model.wSpecies(j)) ')'];
end

%% Constraints?
if isempty(Model.tiepoints); 
    constraints = [];
else
    constraints = ['%% Tiepoints: Depth: ' num2str(Model.tiepoints(1,1)) ...
        'm, Age: ' num2str(Model.tiepoints(1,2)) ' ' Model.ageUnitOut '\r\n'];
    for i = 2:size(Model.tiepoints,1)
        constraints = [constraints '%%            Depth: ' ...
            num2str(Model.tiepoints(i,1)) 'm, Age: ' num2str(Model.tiepoints(i,2)) ' ' Model.ageUnitOut '\r\n'];
    end
end

%% Confidence intervals:
confInterval = [];
nConf = length(Model.confInterval);
for i = 1:nConf
    confInterval = [confInterval num2str(Model.confInterval(nConf-i+1)) ...
        '%% conf (younger bound) \t'];
end
for i = 1:nConf
    confInterval = [confInterval num2str(Model.confInterval(i)) ...
        '%% conf (older bound) \t'];
end
confInterval = confInterval(1:end-2);

%% Header:
header = ['%% Confidence intervals for marker horizons for the ' ...
    Model.icecore ' core \r\n' ...
    '%% ' datestr(now) '\r\n%%\r\n' ...
    '%% Data and weights: ' species_and_weight '\r\n' constraints '%%\r\n'...
    '%% Input layer counts: ' Model.nameManualCounts '\r\n'...
    '%% For method, please cite: \r\n' ...
    '%% Winstrup et al. (2012), Clim. Past 8, 1881-1895, '...
    'An automated method of annual layer detection in ice cores based on '...
    'visual stratigraphy data. \r\n' ...
    '%% Algorithm release date: ' Model.releasedate '\r\n%%\r\n' ...
    '%% Start depth [m] \t End depth [m] \t ML layer \t Prob. of ML layer \t'...
    confInterval '\r\n'];

%% Digits and spacing in data file:
% Number of digits in depth:
if Model.dx_offset == 0
    nDigits = max(ceil(log10(1/Model.dx)));
else
    nDigits = max(ceil(log10(1/Model.dx)))+1;
end

%% Save data with header:
% Open new file with write access, discard any content
fid = fopen([filename '.txt'], 'w');

% Add header:
fprintf(fid, header);
% Add data:
dataformat = ['%.' num2str(nDigits) 'f  \t %.' num2str(nDigits) 'f \t %d \t %.4f'];
for i = 1:nConf
    dataformat = [dataformat '\t %d \t %d'];
end

for iMarkerSet = 1:length(markerConf)
    marker = markerConf{iMarkerSet};
    if ~isempty(marker)
        fprintf(fid, [dataformat '\r\n'], marker');
    end
    % Separate marker sets with header:
    fprintf(fid, '\r\n');
end

% Close file:
fclose(fid);