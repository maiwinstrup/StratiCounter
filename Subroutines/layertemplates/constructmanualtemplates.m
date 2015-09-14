function [Template_all, TemplateInfo_all, Model] = ...
    constructmanualtemplates(Data,Model,outputdir,Runtype)

%% [Template_all, TemplateInfo_all, Model] = ...
%   constructmanualtemplates(Data,Model,outputdir,Runtype)
% Compute array of layer templates for all species corresponding to the 
% manual layer counts in interval "Model.manualtemplates". If running in
% mode Runtype.reuse='yes', layer templates (or preprocessed data to 
% construct these) may be loaded from a previous run. Outputdir is only 
% used when dealing with synthetic data. 

% Copyright (C) 2015  Mai Winstrup
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.

%% Initialize:
switch Model.type
    case 'PCA'
        Template_all(1:Model.nSpecies) = struct('mean',[],'dmean',[],...
            'traj',[],'dtraj',[]);
        TemplateInfo_all(1:Model.nSpecies) = struct('meansignal',[],...
            'pc',[],'score',[],'pcvar',[],'explained',[],'tsquare',[]); 
    case 'FFT'
        Template_all(1:Model.nSpecies) = struct('dc',[],'phase',[],'amplitude',[]);
        TemplateInfo_all(1:Model.nSpecies) = struct('meansignal',[],'meansignal_alt',[]);
end

%% Load manual layer counts for selected depth interval:
[manualcounts, meanLambda, newinterval] = ...
    loadlayercounts(Model,Model.manualtemplates); 

% Changes to depth interval?
if newinterval(1)>Model.manualtemplates(1)
    warning(['Change of interval for calculating layer templates. '...
        'Start: ' num2str(newinterval(1)) 'm'])
end
if newinterval(2)<Model.manualtemplates(2)
    warning(['Change of interval for calculating layer templates. '... 
        'End: ' num2str(newinterval(2)) 'm'])
end
Model.manualtemplates = newinterval;

%% Perform data preprocessing steps for distances given in layer thickness 
% fractions: 
% It is here assumed that the layer thicknesses in the entire interval stay
% relatively constant, so that we can use a fixed processing distance for 
% the data throughout the interval. 
preprocstepsFixed = Model.preprocsteps(:,1);
preprocstepsFloat = setpreprocdist(Model.preprocsteps(:,2),meanLambda);
% Combined data preprocessing steps:
preprocstepsTotal = cell(Model.nSpecies,1);
for j=1:Model.nSpecies
    preprocstepsTotal{j} = [preprocstepsFixed{j}; preprocstepsFloat{j}];
end

%% Load/compute layer templates for the various data records:
for j = 1:Model.nSpecies

    %% Output folder and filenames for layer templates:
    if strcmp(Model.icecore,'SyntheticData')
        filename = [outputdir '/Template_species#' num2str(j)];
    else
        % Template folder for species j:
        % Name describing the preprocessing steps:
        preprocname = makepreprocname(preprocstepsTotal{j},Model.dx);
        % Running in development mode?
        if strcmp(Runtype.develop,'yes'); outputdir0 = './Output/develop';
        else outputdir0 = './Output';
        end        
        
        % Output folder for manual layer templates and parameters:
        outputdir = [outputdir0 '/' Model.icecore '/LayerCharacteristics/' ...
            Model.species{j} '/' preprocname '/' Model.type];
        if ~strcmp(Model.normalizelayer,'none')
            outputdir = [outputdir '_' Model.normalizelayer];
        end
        % Make directory (if doesn't exist):
        if ~exist(outputdir,'dir'); mkdir(outputdir); end

        % Corresponding filename:
        filename = [outputdir '/Template_' num2str(Model.manualtemplates(1)) '-' ...
            num2str(Model.manualtemplates(2)) 'm_' Model.nameManualCounts(1:end-4)];
        if Model.dx_offset~=0
            filename = [filename '_dxoffset' num2str(Model.dx_offset)];
        end
    end

    %% Do layer templates exist from previous run? 
    % If so: Load a previously computed version of templates, unless 
    % a) we're dealing with synthetic data, in which case they are always 
    % recalculated, or b) if Runtype.reuse='no', or c) the file containing 
    % the manual layer counts has been modified since calculation of layer 
    % templates.
    
    % Dates for last modification of files: 
    if exist([filename '.mat'],'file')
        InfoManCounts = dir(['./Manualcounts/' Model.nameManualCounts]);
        dateManCounts = InfoManCounts.datenum; 
        InfoTemplates = dir([filename '.mat']);
        dateTemplates = InfoTemplates.datenum;
    else
        dateManCounts = nan;
        dateTemplates = nan;
    end
    
    if exist([filename '.mat'],'file') && ~strcmp(Model.icecore,'SyntheticData') ...
            && strcmp(Runtype.reuse,'yes') && dateManCounts < dateTemplates
        load([filename '.mat']);
        
        % Check that this template is using the correct value of dtstack:
        if Model.dtstack==1/length(TemplateInfo.meansignal)
            % Use mean signal and principal components from TemplateInfo 
            % to (re)calculate layer templates using the appropriate values 
            % of model order and pcPolOrder. 
            if strcmp(Model.type,'PCA')
                Template_all(j) = polyapprox(TemplateInfo.meansignal,...
                    TemplateInfo.pc,Model); 
            elseif strcmp(Model.type,'FFT')
                disp('FFT version not yet implemented')
            end
            TemplateInfo_all(j) = TemplateInfo;
            disp([Model.species{j} ': Loading old version of layer template'])
            % Proceed to next data record.
            continue
        else
            clear TemplateInfo
            % Recalculate templates corresponding to correct value of
            % dtstack.
        end
    end
       
    %% Otherwise: Calculation of layer shapes!

    %% Preprocess data for interval:
    % Model array corresponding to species j:
    ModelTemplate = Model;
    ModelTemplate.species = Model.species(j); 
    ModelTemplate.nSpecies = 1;
    ModelTemplate.dstart = Model.manualtemplates(1);
    ModelTemplate.dend = Model.manualtemplates(2);
       
    % Does "Data" contain data for the appropriate interval?
    if Data.depth(1)<=Model.manualtemplates(1) && ...
            Data.depth(end)>=Model.manualtemplates(2)
        mask = Data.depth>=Model.manualtemplates(1) & ...
            Data.depth<=Model.manualtemplates(2);
        DataPreproc.depth = Data.depth(mask);
        
        % Finalize preprocessing (if necessary):
        % Using already preprocessed data (which is given as input file):
        DataPreproc.data = makedatafile(Data.data(:,:,j),Data.depth,...
            preprocstepsFloat(j),Model.derivatives,DataPreproc.depth);
        % Derivatives are calculated, but no further downsampling is
        % required.
        
    else
        % Load or, maybe, load and process data for interval.
        % Processed data may not exist from before if e.g. the provided 
        % layer template interval is not covered by data interval for layer
        % counting. 
        
        % Construct processed data file: 
        % All preprocessing steps must be performed. 
        ModelTemplate.preprocsteps = preprocstepsTotal(j);
        DataPreproc = loadormakedatafile(ModelTemplate,manualcounts,Runtype); 
    end
    
    %% Calculate layer templates:
    [Template, TemplateInfo] = layerstructure(DataPreproc.data,...
        DataPreproc.depth,manualcounts(:,1),manualcounts(:,3),...
        ModelTemplate,Runtype); 
    
    %% Save templates:
    % As data file:
    save([filename '.mat'],'TemplateInfo')
    
    % ...and in array "Template_all":
    Template_all(j)=Template;
    TemplateInfo_all(j)=TemplateInfo;
    
    %% Plot principal components, and save figure:
    if Runtype.plotlevel > 0
        if strcmp(Model.type,'PCA')
            plotprincomp(TemplateInfo,Model.species{j},Model.icecore,Model.nameManualCounts)
            print([filename '.jpeg'], '-djpeg','-r400')
            % Close plot (if plotlevel<2):
            if Runtype.plotlevel==1; close(gcf); end
        end
    end
end
end

function plotprincomp(TemplateInfo,species,icecore,manCountsName)
% Plot mean signal and principal components of layers of "species".
% Mai Winstrup, 2015

%% Mean signal:
figure;
subplot(2,1,1)  
x = 1/(2*length(TemplateInfo.meansignal)):1/length(TemplateInfo.meansignal):1;
plot(x,TemplateInfo.meansignal,'-k','linewidth',2)
% Add label and title:
ylabel('Mean signal','fontweight','bold')
if strcmp(icecore,'SyntheticData')
    title(['Species #' species ' (' manCountsName ')'],'fontweight','bold','interpreter','none')
else
    title([species ' (' manCountsName ')'],'fontweight','bold','interpreter','none')
end

%% Principal components (plotting the first three components): 
subplot(2,1,2)
plot(x,TemplateInfo.pc(:,1:3),'linewidth',2)
ylabel('PCs','fontweight','bold')
% Add legend:
for i = 1:3; legendname{i,:} = ['PC' num2str(i)]; end
legend(legendname)
legend('boxoff')
end