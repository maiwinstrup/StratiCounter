function sectionswithoutdata(Data,sectionlength,species)

%% Check for long sections without data:
% "Long sections without data" is defined as sections larger than 
% "sectionlength" for which less than half of section contains data.
% Copyright (C) 2015  Mai Winstrup

%% Section length in pixels:
dx = mean(diff(Data.depth)); % Using an average value for dx
if isnan(sectionlength)
    sectionlength = input('Length of section (in m) defining a large data gab (usually 20*lambda): ');
end
sectionlength_px = round(sectionlength/dx); %[pixel]

%% Find nans in data (not considering derivatives) 
% And filter the resulting logical array.
nSpecies = size(Data.data,3);
nanData = nan(length(Data.depth),nSpecies);
for j = 1:nSpecies
    mask = isnan(Data.data(:,1,j)); % not considering derivatives
    nanData(:,j)=medfilt1(mask*1,sectionlength_px);
end

% Number of missing data series for each pixel:
nNanData = sum(nanData,2);

%% Start and end boundaries of these sections:
% Defined as areas where the number of nan data series increase/decrease, 
% *and* start from the baseline of zero (i.e. all data exist). 
nanDataBounds = diff(nNanData);
istart = find(nanDataBounds>0);
% And starting from zero:
mask = nNanData(istart)==0;
istart = istart(mask);
% Inserting an extra boundary in the beginning if starting out with nans:
if nNanData(1)~=0; istart = [1; istart]; end
dstart = Data.depth(istart);

% Ending boundaries:
iend = find(nanDataBounds<0);
% And ending with no data series with nan:
mask = nNanData(iend+1)==0;
iend = iend(mask);
% Inserting an extra boundary if ending with nans:
if nNanData(end)~=0; iend = [iend; length(Data.depth)]; end
dend = Data.depth(iend);

% Length of sections:
Lnan = dend-dstart;

% Number of data series with nan in each section - and which ones:
for i = 1:length(dstart)
    % Maximum number of data series with nan:
    nNanSpecies(i) = max(nNanData(istart(i):iend(i)));
    % The data series are:
    mask = sum(nanData(istart(i):iend(i),:),1)>0;
    nanSpecies{i} = species(mask);
end

% Display error message, sorted relative to depth: 
if sum(nNanData>0)
    disp('OBS: The following are large sections without much data:')
    for i = 1:length(dstart)
        % Create a list of missing data files in section:
        species = [];
        for k = 1:length(nanSpecies{i})
            species = [species nanSpecies{i}{k} ', '];
        end
        species = species(1:end-2); % removing last comma
        % Display message:
        if nNanSpecies(i)==nSpecies
            disp([num2str(dstart(i)) '-' num2str(dend(i)) '(' num2str(Lnan(i))...
            'm): All data series missing'])
        else
            disp([num2str(dstart(i)) '-' num2str(dend(i)) '(' ...
                num2str(Lnan(i)) 'm): ' num2str(nNanSpecies(i)) ' data '...
                'series missing (' species ')'])
        end
    end 
end