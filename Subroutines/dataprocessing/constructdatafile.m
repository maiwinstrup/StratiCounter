function [Data_final, Model] = constructdatafile(Model,counts,Runtype)

%% [Data_final, Model] = constructdatafile(Model,counts,Runtype)
% Load and/or process selected impurity records for ice core. Preprocessing 
% is done according to Model.preprocessing_init. Processed data is saved. 
% If Runtype.reuse='no', preprocessing is done from scratch, otherwise 
% processed data may be loaded from prior runs. Counts are only used for 
% plotting, and do not need to be specified. 

% Copyright (C) 2015  Mai Winstrup
% 2014-04-21 17:15: Filename updated, figures may close after use
% 2014-06-24 14:17: Runtype.startover -> Runtype.reuse
% 2014-07-15 10:09: New structure of Model.preprocess, and use of floating distances
%                   Changed depth scale for data series with "half-way" points
% 2014-07-16 16:50: Fixed bug reg. double-processing of non-equidistant
%                   data, non-equidistant data allowed
% 2014-07-24 11:06: Solved problem with constant reprocessing of data series, which
%                   was caused by depth scale inaccuracies.
% 2014-08-06 17:16: Solved problem with depth scale inaccuracies,
%                   dx_midpoint -> dx_center
% 2014-08-15 17:13: Removal of Runtype.test (using develop instead)
% 2014-08-16 13:37: Name changed of createpreprocessname
% 2014-08-17 22:01: fixed "minusmean" i processdata.
% 2014-08-21 22:56: function to construct outputdir introduced, changes to
%                   downsampling, ensure correct values of dstart and dend
% 2014-08-22 14:52: correct values of dstart and dend moved to checkmodel
% 2014-08-23 13:58: few changes to processing of data 
% 2014-09-08 14:02: removed +dx_center in line 72
% 2014-10-01 13:06: solved problem when not all data series cover the
%                   entire interval
% 2015-01-21 12:17: Changes due to new definition of dx_center, new
%                   datafile structure

%% Provide indicator of missing impurity records:
flag = nan(1,Model.nSpecies);

% Final depth scale for data:
Data_final.depth(:,1) = Model.dstart+Model.dx_center*Model.dx:Model.dx:Model.dend+Model.dx_center*Model.dx;
Data_final.data = nan(length(Data_final.depth),3,Model.nSpecies);

%% Load/process the individual data records:
for j = 1:Model.nSpecies

    %% Output folder for processed data:
    outputdir = makepreprocfolder(Model.preprocess_init{j},Model.dx,Model.icecore,...
        Model.species{j},Runtype); % 2014-08-21 21:12
    
    %% Does processed data exist? 
    % If so, these are loaded, provided Runtype.reuse='yes'
    if strcmp(Runtype.reuse,'yes') 
        % List existing datafiles with similar resolution and preprocessing:
        listing = dir([outputdir '/*m.mat']);

        % Removing those not covering the appropriate depth interval:
        N = size(listing,1);
        if N>0
            d1 = zeros(1,N);
            d2 = zeros(1,N);
            for i = 1:N
                index1 = strfind(listing(i,1).name,'-');
                index2 = strfind(listing(i,1).name,'m');
                d1(i) = str2double(listing(i,1).name(1:index1-1));
                d2(i) = str2double(listing(i,1).name(index1+1:index2-1));
            end
            mask = d1<=Model.dstart & d2>=Model.dend-Model.dx/100;
            listing = listing(mask,:);
        
            % Load the smallest file (if file exist):
            if size(listing,1)>0
                [~,imin]=min([listing.bytes]);
                filename = listing(imin,1).name;
                load([outputdir '/' filename]);
        
                % Remove excess data from outside interval:
                mask = depth>=Model.dstart+Model.dx_center*Model.dx & depth<=Model.dend+Model.dx_center*Model.dx;
                
                % Check that depth scale is correct, which can be seen from 
                % the initial depth entry. If not, the data will be 
                % re-processed to correct depth scale.
                depth = depth(mask);
                if depth(1)==Model.dstart+Model.dx_center*Model.dx
                    % Add to output data file:
                    Data_final.data(:,:,j) = data(mask,:);
                
                    % Test for only NaNs in data series:
                    if sum(isfinite(Data_final.data))==0; 
                        flag(j) = 1; % Data are nan in all of current interval 
                        disp([model.species{j} ' not available in interval'])
                    end
                
                    % Also loaded are the analytically-derived relative 
                    % white noise weighting values for the derivative data 
                    % series. 
                    % If appropriate, these are added to the Model structure. 
                    if strcmp(Model.wWhiteNoise,'analytical')
                       Model.wWhiteNoise = wWhiteNoise;
                    end
                
                    % And we're done! 
                    % Proceed to next impurity record.
                    continue
                else
                    clear data depth
                end
            end
        end
    end

    %% Otherwise, load and preprocess the raw data for current interval.
    rawdata = loadrawdata(Model.species{j},Model); % 2015-01-21

    %% If data record does not exist:
    % Set flag and go to next impurity species.
    if isempty(rawdata); 
        flag(j) = 2; disp([Model.species{j} ' do not exist']);
        continue
    end
    
    % Does data exist in at least part of data interval?   
    % If not, set flag, and continue to next impurity species.
    mask = rawdata(:,1)>=Model.dstart+Model.dx_center*Model.dx & rawdata(:,1)<=Model.dend+Model.dx_center*Model.dx;
    if sum(isfinite(rawdata(mask,2)))==0; 
        % Data is nan in all of current interval 
        flag(j) = 1; disp([Model.species{j} ' not available in interval']); 
        continue 
    end
        
    %% Ensure that the data series is sorted according to depth, and that 
    % only a single data point corresponds to a specific depth:
    rawdata = sortrows(rawdata,1);

    % Some data sets have multiple data values for a single depth entry. 
    % These are removed and replaced by an average value:
    index=find(diff(rawdata(:,1)==0));
    if ~isempty(index)
        for i = index(1):index(end)
            % There may be multiple datapoints in a row assigned the same 
            % depth:
            mask = rawdata(:,1)==rawdata(i,1);
            % All are given the mean value
            meanvalue = mean(rawdata(mask,2));
            rawdata(mask,2) = meanvalue;
        end
        % Remove double entries:
        rawdata(index,:)=[];
    end
    
    %% Remove data from outside interval, while keeping some extra data 
    % around edges. 
    mask = rawdata(:,1)>=Model.dstart-10 & rawdata(:,1)<=Model.dend+10;
    depth = rawdata(mask,1);
    data0 = rawdata(mask,2);
    
    %% Preprocess data:
    [data0,hfig(1),dstart_fig,dend_fig] = processdata(data0,Model.species{j},...
        Model.preprocess_init{j},Model,depth,counts,Runtype.plotlevel); % 2014-08-23 13:55
    
    %% Remove extra data from edges: 
    mask = depth>=Model.dstart+Model.dx_center*Model.dx & depth<=Model.dend+Model.dx_center*Model.dx;
    depth = depth(mask);
    data0 = data0(mask,:);
    
    %% Downsample to resolution in Model.dx:
    % If Model.dx is empty, the original resolution is kept.
    if ~isempty(Model.dx)
        [depth,data0] = downsampling(depth,data0,Model.dx,Model.dx_center*Model.dx,Model.dstart,Model.dend); % % 2014-08-21 22:33
    end
    
    %% Calculating derivatives of data series:
    [slope, dslope, wWhiteNoise] = ...
        calculateslope(data0,Model.slopeorder,Model.slopedist,0); % 2014-07-16 14:05
    % The derivatives are calculated per pixel, and thus their calculation 
    % do not require an equidistant timescale.

    % Make data file:
    clear data
    data(:,1) = data0;
    data(:,2) = slope;
    data(:,3) = dslope;
    
    % Plot derivatives:
    if Runtype.plotlevel>0
        hfig(2) = plotderivatives(depth,data,dstart_fig,dend_fig,Model,j,counts);
    end

    %% Save data and figures:
    filename = [num2str(Model.dstart) '-' num2str(Model.dend) 'm'];
    save([outputdir '/' filename '.mat'],'data','depth','wWhiteNoise');
    
    % Save example of datafiles:
    if Runtype.plotlevel>0
        figure(hfig(1));
        print([outputdir '/' filename '_data.jpeg'],'-djpeg','-r300')
        figure(hfig(2));
        print([outputdir '/' filename '_deriv.jpeg'],'-djpeg','-r300')
    end

    % And (maybe) close figures:
    if Runtype.plotlevel==1
        close(hfig(1)); close(hfig(2))
    end
    
    %% Include in data array containing data for all species:
    % Place data correctly in array:
    index = interp1(Data_final.depth,1:length(Data_final.depth),depth,'nearest');
    Data_final.data(index,:,j) = data; 
    
    % Display warning if data does not cover all of interval:
    if index(1)>10 || index(end) < length(Data_final.depth)-10
        disp(['OBS: ' Model.species{j} ' only covers interval ' num2str(depth(1)-Model.dx_center) '-' num2str(depth(end)-Model.dx_center) 'm.'])
    end
    
    %% Add analytical weighting values to Model (if wanted):
    if strcmp(Model.wWhiteNoise,'analytical')
        Model.wWhiteNoise = wWhiteNoise;
    end
end

%% Remove unavailable data records from consideration as data series for 
% the annual layering:
mask = isfinite(flag);
if sum(mask)>0
    Data_final.data(:,:,mask) = [];
    Model.species = Model.species(~mask);
    Model.nSpecies = length(Model.species);
    Model.preprocess = Model.preprocess(~mask);
    Model.preprocess_init = Model.preprocess_init(~mask);
    Model.preprocess_rep = Model.preprocess_rep(~mask);
    Model.wSpecies = Model.wSpecies(~mask);
end

%% Check for long sections without data:
% "Long sections without data" is defined as sections larger than ~N years 
% for which less than half of section contains data.
N = 20; % Layers
% Length of sections: 
meanlambda = mean(diff(counts(:,1)));
L = N*meanlambda; % [m]
Lpx = round(L/Model.dx); %[pixel]

nanData = nan(length(Data_final.depth),Model.nSpecies);
for j = 1:Model.nSpecies
    mask = isnan(Data_final.data(:,1,j)); % not looking at derivatives
    nanData(:,j)=medfilt1(mask*1,Lpx);
end
% Number of missing data series for each pixel:
nNanData = sum(nanData,2);

% Start and end boundaries of these sections:
% Defined as areas where the number of nan data series increase/decrease, 
% *and* start from the baseline of zero (i.e. all data exist). 
nanDataBounds = diff(nNanData);
istart = find(nanDataBounds>0);
% And starting from zero:
mask = nNanData(istart)==0;
istart = istart(mask);
% Inserting an extra boundary in the beginning if starting out with nans:
if nNanData(1)~=0; istart = [1; istart]; end
dstart = Data_final.depth(istart);

% Ending boundaries:
iend = find(nanDataBounds<0);
% And ending with no data series with nan:
mask = nNanData(iend+1)==0;
iend = iend(mask);
% Inserting an extra boundary if ending with nans:
if nNanData(end)~=0; iend = [iend; length(Data_final.depth)]; end
dend = Data_final.depth(iend);

% Length of sections:
Lnan = dend-dstart;

% Number of data series with nan in each section - and which ones:
for i = 1:length(dstart)
    % Maximum number of data series with nan:
    nNanSpecies(i) = max(nNanData(istart(i):iend(i)));
    % The data series are:
    mask = sum(nanData(istart(i):iend(i),:),1)>0;
    nanSpecies{i} = Model.species(mask);
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
        if nNanSpecies(i)==Model.nSpecies
            disp([num2str(dstart(i)) '-' num2str(dend(i)) '(' num2str(Lnan(i))...
            'm): All data series missing'])
        else
            disp([num2str(dstart(i)) '-' num2str(dend(i)) '(' ...
                num2str(Lnan(i)) 'm): ' num2str(nNanSpecies(i)) ' data '...
                'series missing (' species ')'])
        end
    end 
end
end

%% Plotting subroutine:
function hfig = plotderivatives(depth,data_final,dstart_fig,dend_fig,Model,j,counts)
hfig = figure;
for i = 1:3
    subplot(3,1,i); 
    plot(depth,data_final(:,i))
    xlim([dstart_fig dend_fig])
    hold on
    % Plot layer positions:
    plotlayercounts(counts,data_final(:,i))
end
subplot(3,1,1)
title(['Data: ' Model.species{j}],'fontweight','bold'); 
subplot(3,1,2) 
title(['Slope (order=' num2str(Model.slopeorder(1)) ', L=' num2str(Model.slopedist(1)) ')'],'fontweight','bold')
subplot(3,1,3)
title(['Curvature (order=' num2str(Model.slopeorder(2)) ', L=' num2str(Model.slopedist(2)) ')'],'fontweight','bold')   
end