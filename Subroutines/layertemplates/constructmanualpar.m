function Layerpar0 = constructmanualpar(Data,Template,Model,outputdir,Runtype)

%% Layerpar0 = constructmanualpar(Data,Template,Model,outputdir,Runtype)
% Calculate layer parameters corresponding to the given template, and
% manual layer counts. Layer parameters are always recalculated to prevent 
% troubles with matching to adjusted templates.
% Copyright (C) 2015  Mai Winstrup
 
%% Load manual layer counts for selected depth interval:
[manualcounts,meanLambda,newinterval] = ...
    loadlayercounts(Model,Model.initialpar);
if isempty(manualcounts)
    error(['Manual counts do not exist in selected interval for initial '...
        'parameter estimation.'])
end

% Changes to depth interval?
if newinterval(1)>Model.initialpar(1)
    warning(['Note: Change of interval for calculating initial layer parameters. '...
        'Start: ' num2str(newinterval(1)) 'm'])
end
if newinterval(2)<Model.initialpar(2)
    warning(['Note: Change of interval for calculating initial layer parameters. '... 
        'End: ' num2str(newinterval(2)) 'm'])
end
Model.initialpar = newinterval;

%% Perform data preprocessing steps for distances given in layer thickness 
% fractions. If depth interval is outside data range, data for interval is 
% first loaded and then preprocessed according to all preprocessing steps. 

% For the floating preprocessing steps it is assumed that the layer 
% thicknesses in the entire interval stay relatively constant, so that we 
% can use a fixed processing distance for the data throughout the interval. 
preprocstepsFloat = setpreprocdist(Model.preprocsteps(:,2),meanLambda);

% If initial parameter interval is within depth section for layer counting: 
if Model.initialpar(1)>=Data.depth(1) && Model.initialpar(2)<=Data.depth(end)
    % Data in interval:
    mask=Data.depth>=Model.initialpar(1)&Data.depth<=Model.initialpar(2);
    depth = Data.depth(mask);
    data = Data.data(mask,:,:);
    
else
    % Load and/or process data for interval.
    % Processed data may not exist from before if e.g. the provided 
    % parameter interval is not covered by data interval for layer counting.  
        
    % Construct processed data file: 
    % No plotting:
    Runtype1 = Runtype; Runtype1.plotlevel=0;
    % For layerpar interval:
    Model1 = Model; 
    Model1.dstart = Model.initialpar(1);
    Model1.dend = Model.initialpar(2);
    % Load and/or process data:
    Data = loadormakedatafile(Model1,[],Runtype1);
    depth = Data.depth;
    data = Data.data;
end
    
% Finalize the preprocessing steps:
DataPreproc.depth = depth;
DataPreproc.data = makedatafile(data,depth,preprocstepsFloat,...
    Model.derivatives,depth); 
% Keeping the same depth scale - no further downsampling etc.

%% Compute layer parameters for the various data records:
% Calculating both the maximum-Likelihood initial parameters (i.e. the best 
% fitting parameters for each layer, based on both data and derivatives)
% and the Maximum-a-Posteriori layer parameters: 
Layerpar0 = calclayerpar(Model,DataPreproc,manualcounts(:,1),...
    manualcounts(:,3),Template,Runtype);

%% Not using full covariance matrix?
if strcmp(Model.covariance,'none')
    Layerpar0.cov = diag(diag(Layerpar0.cov));
end

%% If no parameter updates: 
% Ask user if desire to use the values from manual counts.
noupdates = strcmp(Model.update,'none');
index = find(noupdates==1);
names = {'my','sigma','par','cov','nvar'};
for i= 1:length(index)
    disp(['The value of ' names{index(i)} ' will be held constant.'])
    disp('Value corresponding to manual counts is: ');
    disp(num2str(eval(['Layerpar0.' names{index(i)}])))
    reply = input('Use this value? (y/n) ','s');
    if strcmp(reply,''); disp('yes'); end
    
    if ~ismember(reply,{'yes','y','Y','YES',''})       
        value = input('Select new value: ');
        % Check for correct format:
        while ~isequal(size(value),size(eval(['Layerpar0.' names{index(i)}])))
            value = input('Incorrect format, please correct: ');
        end
        while (ismember(names{index(i)},{'sigma','nvar'})&& sum(value<=0)>0) || ...
                strcmp(names{index(i)},{'cov'}) && sum(diag(value)<=0)>0
            value = input(['Value of ' names{index(i)} ' must be great than 0. Please correct: ']);
        end
        % Inset value in array:
        if index(i) == 1; Layerpar0.my = value;
        elseif index(i) == 2; Layerpar0.sigma = value;
        elseif index(i) == 3; Layerpar0.par = value;
        elseif index(i) == 4; Layerpar0.cov = value;
        elseif index(i) == 5; Layerpar0.nvar = value;
        end
        disp('New value is: '); disp(value)
    end
end

%% Save the employed set of layerparameters in timescale output folder:
save([outputdir '/Layerpar0'],'Layerpar0')