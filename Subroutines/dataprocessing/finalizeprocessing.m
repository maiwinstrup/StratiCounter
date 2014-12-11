function [depth_interval,data_final] = finalizeprocessing(depth,data,meanLambda,interval,Model,loglevel)

%% [depth_interval, data_final] = finalizeprocessing(depth,data,meanLamda,
% interval,Model,loglevel)
% Perform the final processing according to the floating distances, if such
% ones exist. Format of data_final is the same as for data.

% Mai Winstrup
% 2014-08-23 14:50: Initial version
% 2014-10-01 13:19: 

%% Data in interval:
mask=depth>=interval(1)&depth<=interval(2);
depth_interval = depth(mask);
data_interval = data(mask,:,:);

%% Process data:
% Data need further processing if we have floating distances. 
% Any remaining processing steps are contained in Model.preprocess_rep, and
% the corresponding numeric values of the preprocessing distances are:
preprocess_rep = setfloatingdist(Model.preprocess_rep,meanLambda); % 2014-08-16 12:27

% Process data, if necessary:
data_final = nan(size(data_interval));
for j = 1:Model.nSpecies
    % Check if data record is nan in entire interval:
    finitedata = sum(isfinite(data(:,1,j))); 
    if finitedata(1)==0
        % No need for further processing (it's all NaN anyway):
        data_final(:,:,j) = data_interval(:,:,j);
        % Display error message: 
        if loglevel>0
            disp([Model.species{j} ' not available in interval: ' num2str(interval) 'm'])
        end
        
    elseif strcmp(preprocess_rep{j}(1,1),'none')
        % No floating distances, no further processing: 
        data_final(:,:,j) = data_interval(:,:,j);
        
    else
        % Using a slightly extended area: 
        % Additional number of datapoints:
        Lext = max(cell2mat(preprocess_rep{j}(:,2)))/Model.dx;
        if isempty(Lext); Lext = 0; end
        
        % Start and end of data section to be used: forkert!
        i0start = find(depth==depth_interval(1),1,'first');
        istart = max(1,i0start);
        i0end = find(depth==depth_interval(end),1,'first');
        iend = min(length(depth),i0end);
        
        % Process data in extended area:
        data0 = processdata(data(istart:iend,1,j),Model.species{j},...
            preprocess_rep{j},Model,depth(istart:iend)); % 2014-08-23 13:55
        % Remove extended area:
        data0 = data0(1+(i0start-istart):end-(iend-i0end)); 
        
        % Add to resulting data file:
        data_final(:,1,j) = data0;

        if size(data,2) > 1
            % Calculate derivatives: 
            [slope, dslope] = calculateslope(data0,Model.slopeorder,...
                Model.slopedist,0); % 2014-07-16 14:05
            data_final(:,2,j) = slope;
            data_final(:,3,j) = dslope;
        end
    end
end