function rawdata = loadrawdata(species,path2data)

%% data = loadrawdata(species,path2data)
% Load the corresponding impurity record.
% Output:   rawdata(:,1): Depth
%           rawdata(:,2): Data
% The data is sorted depthwise, and checked for multiple data values 
% corresponding to a single depth entry. 
% Copyright (C) 2015  Mai Winstrup

%% Default output:
rawdata = [];

%% Load data for ice core:
if ~exist(path2data,'file')
    disp('Path to datafile is not working')

else
    % Load data structure array:
    load(path2data)
    
    % Find location of species in array: 
    index = find(strcmp(data.name,species)==1);
    if isempty(index)
        disp([species ' does not exist in data file'])
        return
    end
        
    % Corresponding depth scale:
    depth_no = data.depth_no(index);
    rawdata(:,1) = data.depth{depth_no};

    % Data file:
    rawdata(:,2) = data.data{index};
end

%% Sort data record according to depth:
rawdata = sortrows(rawdata,1);

%% Ensure that only a single data point corresponds to a specific depth:
% Some data sets may have multiple data values for a single depth entry. 
% These are removed and replaced by an average value:
index=find(diff(rawdata(:,1)==0));
if ~isempty(index)
    for i = index(1):index(end)
        % There may be multiple datapoints in a row assigned the same depth:
        mask = rawdata(:,1)==rawdata(i,1);
        % All are given the mean value:
        meanvalue = mean(rawdata(mask,2));
        rawdata(mask,2) = meanvalue;
    end
    % Remove double entries:
    rawdata(index,:)=[];
end