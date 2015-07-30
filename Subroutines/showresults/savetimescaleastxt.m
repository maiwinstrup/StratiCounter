function savetimescaleastxt(timescale1yr,filename,Model)

%% savetimescaleastxt(timescale1yr,filename,Model)
% Save timescale results as textfile with metadata and header.
% Copyright (C) 2015  Mai Winstrup

%% Meta data:
% Species and their associated weighting:
species_and_weight = [Model.species{1} ' (w=' num2str(Model.wSpecies(1)) ')'];
for j = 2:Model.nSpecies
    species_and_weight = [species_and_weight ', ' Model.species{j}  ...
        ' (w=' num2str(Model.wSpecies(j)) ')'];
end

% Constraints?
if isempty(Model.tiepoints); 
    constraints = [];
else
    constraints = ['%% Tiepoints: Depth: ' num2str(Model.tiepoints(1,1)) ...
        'm, Age: ' num2str(Model.tiepoints(1,2)) ' ' Model.ageUnitOut '\r\n'];
    for i = 2:size(Model.tiepoints,1)
        constraints = [constraints '%%            Depth: ' ...
            num2str(Model.tiepoints(i,1)) 'm, Age: ' num2str(Model.tiepoints(i,2)) ...
            ' ' Model.ageUnitOut '\r\n'];
    end
end

% Confidence intervals:
confInterval = [];
nConf = length(Model.confInterval);
for i = 1:nConf
    confInterval = [confInterval num2str(Model.confInterval(nConf-i+1)) ...
        '%% conf int (younger bound) \t '];
end
for i = 1:nConf
    confInterval = [confInterval num2str(Model.confInterval(i)) ...
        '%% conf int (older bound) \t '];
end
confInterval = confInterval(1:end-3);

%% Construct metadata and header:
header = ['%% Timescale for the ' Model.icecore ' core \r\n' ...
    '%% ' datestr(now) '\r\n%%\r\n' ...
    '%% Data and weights: ' species_and_weight '\r\n' constraints '%%\r\n'...
    '%% Input layer counts: ' Model.nameManualCounts '%%\r\n'...
    '%% For method, please cite: \r\n' ...
    '%% Winstrup et al. (2012), Clim. Past 8, 1881-1895, '...
    'An automated method of annual layer detection in ice cores based on '...
    'visual stratigraphy data. \r\n' ...
    '%% Algorithm release date: ' Model.releasedate '\r\n%%\r\n' ...
    '%% Depth [m] \t Maximum Likelihood age [' Model.ageUnitOut '] \t '...
    confInterval ' \r\n'];

%% Digits and spacing in data file:
% Number of digits in depth:
if Model.dx_offset == 0
    nDigits = max(ceil(log10(1/Model.dx)));
else
    nDigits = max(ceil(log10(1/Model.dx)))+1;
end

%% Save data with header:
% Open new file with write access, discard any content
fid = fopen(filename, 'w');

% Add header:
fprintf(fid, header);
% Add data:
timescaleformat = ['% .' num2str(nDigits) 'f  \t %d'];
for i = 1:nConf
    timescaleformat = [timescaleformat '\t %d \t %d'];
end
fprintf(fid, [timescaleformat '\r\n'], timescale1yr'); % First: space, second: #digits

% Close file:
fclose(fid);