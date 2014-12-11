function rawdata = loadrawdata(species,Model)

%% data = loadrawdata(species,Model)
% Load corresponding impurity record.
% Output:   data(:,1): Depth
%           data(:,2): Data
% Mai Winstrup
% 2014-04-21 17:03: Filename updated

%% Set as default:
rawdata = [];

%% Load data for ice core:
if ~exist(Model.pathData,'file')
    disp('Path to datafile is not working')

else
    load(Model.pathData)

    % Depth scale:
    rawdata(:,1) = data.depth;

    % Data file:
    index = find(strcmp(data.name,species)==1);
    if isempty(index)
        disp([species ' does not exist in data file'])
    else
        rawdata(:,2) = data.data(:,index);
    end
end