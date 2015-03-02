function Layerpar0 = constructmanualpar(Data,Template,Model,outputdir,Runtype)

%% Layerpar0 = constructmanualpar(Data,Template,Model,outputdir,Runtype)
% Calculate layer parameters corresponding to the given template, and
% manual layer counts. Layer parameters are always recalculated to prevent 
% troubles with matching to adjusted templates.
% Copyright (C) 2015  Mai Winstrup
   
% ville det være smart at det ikke nødvnedigvis var baseret på man counts?   
% de gemmes jo alligevel ikke. 
% ønsker vi f.eks at benytte den i loop, hvor vi har udregnet nye counts,
% og starter forfra, og vil gerne have et par gode parametre som passer til
% vores nye template?
 
% i virkeligheden: 2 fsk. een for manuelt,  som kalder en funktion som
% udregner det baseret på data og "counts", for et givet interval. Denne
% kan bruges senere. 


%% Load manual layer counts for selected depth interval:
[manualcounts,meanLambda,newinterval] = ...
    loadlayercounts(Model,Model.initialpar); % 2014-08-22 16:30

% Changes to depth interval?
if newinterval(1)>Model.initialpar(1)
    disp(['Note: Change of interval for calculating initial layer parameters. '...
        'Start: ' num2str(newinterval(1)) 'm'])
end
if newinterval(2)<Model.initialpar(2)
    disp(['Note: Change of interval for calculating initial layer parameters. '... 
        'End: ' num2str(newinterval(2)) 'm'])
end
Model.initialpar = newinterval;

%% Processed data for interval:
% Initial parameter interval is assumed to always be within the depth range 
% of our depth-interval, i.e. contained in Data. 
% Data in interval:
mask=Data.depth>=Model.initialpar(1)&Data.depth<=Model.initialpar(2);
depth = Data.depth(mask);
data = Data.data(mask,:,:);
       
% Finalize the processing:
preprocsteps=setpreprocdist(Model.preprocess(:,2),meanLambda);
[DataProcessed.data, DataProcessed.depth] = makedatafile(data,depth,...
    preprocsteps,Model.derivatives); % No further downsampling or plotting

%% Compute layer parameters for the various data records:
% Calculating both the maximum-Likelihood initial parameters (i.e. the best 
% fitting parameters for each layer, based on both data and derivatives),
% and the Maximum-a-Posteriori layer parameters: 
[~,parML,parMAP,Layerpar0] = calclayerpar(Model,DataProcessed,manualcounts(:,1),...
    manualcounts(:,3),Template,Runtype);

%% Not using full covariance matrix?
if strcmp(Model.covariance,'none')
    Layerpar0.cov = diag(diag(Layerpar0.cov));
end

%% Save the employed set of layerparameters in output folder for timescale:
save([outputdir '/Layerpar0'],'Layerpar0')

%         % OBS: hvad så med Model.w - er ok fra tidligere????