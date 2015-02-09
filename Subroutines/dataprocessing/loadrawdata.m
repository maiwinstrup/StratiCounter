function rawdata = loadrawdata(species,Model)

%% data = loadrawdata(species,Model)
% Load corresponding impurity record.
% Output:   data(:,1): Depth
%           data(:,2): Data

% Copyright (C) 2015  Mai Winstrup
% 2014-04-21 17:03: Filename updated
% 2015-01-21 11:44: New data structure array

%% Set as default:
rawdata = [];

%% Load data for ice core:
if ~exist(Model.pathData,'file')
    disp('Path to datafile is not working')

else
    load(Model.pathData)
    index = find(strcmp(data.name,species)==1);

    if isempty(index)
        disp([species ' does not exist in data file'])
        return
    end
        
    % Depth scale:
    depth_no = data.depth_no(index);
    rawdata(:,1) = data.depth{depth_no};

    % Data file:
    rawdata(:,2) = data.data{index};
end