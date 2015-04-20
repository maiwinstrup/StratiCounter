function [Template_all, TemplateInfo_all] = ...
    constructmanualtemplates(Data,Model,outputdir,Runtype)

%% [Template_all, TemplateInfo_all] = constructmanualtemplates(Data,Model,outputdir,Runtype)
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

% 2014-04-21 23:21: Small adjustments
% 2014-08-15 15:40: New name for function, shape->template
% 2014-08-17 22:46: Changes to layerstructure etc.
% 2014-08-18 15:10: Data as input argument to function to solve problems 
%                   with possible double recalculation when reuse='no',
%                   input arguments to layer structure changed.
% 2014-08-19 23:25: Changes to load of layercounts and outputdir for 
%                   synthetic data, changes to layerstructure
% 2014-08-21 14:11: Introduced outputdir for synthetic data, and using
%                   maketemplate for normal data.
% 2014-08-21 23:22: new structure of loadlayercounts, changes to
%                   constructdatafile
% 2014-08-22 16:30: small changes
% 2014-10-01 13:21: still problems when numeric value in preproc does not 
%                   provide a distance, but some other measure. Issue not solved.

%% Initialize::
switch Model.type
    case 'PCA'
        Template_all(1:Model.nSpecies) = struct('mean',[],'dmean',[],...
            'traj',[],'dtraj',[]);
        TemplateInfo_all(1:Model.nSpecies) = struct('meansignal',[],...
            'pc',[],'score',[],'pcvar',[],'explained',[],'tsquare',[]); 
    case 'FFTcomp'
        Template_all(1:Model.nSpecies) = struct('dc',[],'phase',[],'amplitude',[]);
        TemplateInfo_all(1:Model.nSpecies) = struct('meansignal',[],'meansignal_alt',[]);
end

%% Load manual layer counts for selected depth interval:
[manualcounts, meanLambda, newinterval] = ...
    loadlayercounts(Model,Model.manualtemplates); % 2014-08-23 13:21

% Changes to depth interval?
if newinterval(1)>Model.manualtemplates(1)
    disp(['Note: Change of interval for calculating layer templates. '...
        'Start: ' num2str(newinterval(1)) 'm'])
end
if newinterval(2)<Model.manualtemplates(2)
    disp(['Note: Change of interval for calculating layer templates. '... 
        'End: ' num2str(newinterval(2)) 'm'])
end
Model.manualtemplates = newinterval;

%% Set data pre-processing distances for templates: 
% Convert from floating distances to numeric values. We use a fixed 
% distance for the data series throughout the interval. 
preprocessTemplate = setpreprocdist(Model.preprocess(:,2),meanLambda);

%% Load/compute layer templates for the various data records:
for j = 1:Model.nSpecies

    %% Output folder and filenames for layer templates:
    if strcmp(Model.icecore,'SyntheticData')
        filename = [outputdir '/Template_species#' num2str(j)];
    else
        % Template folder for species j:
        outputdir = maketemplatefolder(preprocessTemplate{j},Model.dx,...
            Model.icecore,Model.species{j},Model.type,...
            Model.normalizelayer,Runtype); % 2014-08-22 15:00
        filename = [outputdir '/Template_' num2str(Model.manualtemplates(1)) '-' ...
            num2str(Model.manualtemplates(2)) 'm'];
    end

    %% Do layer templates exist from previous run? 
    % If so: Load a previously computed version of templates, unless 
    % a) we're dealing with synthetic data, in which case they are always 
    % recalculated, or b) if Runtype.reuse='no'.    
    
    if exist([filename '.mat'],'file') && ~strcmp(Model.icecore,'SyntheticData') ...
            && strcmp(Runtype.reuse,'yes')
        load([filename '.mat']);
        
        % Check that this template is using the correct value of dtstack:
        if Model.dtstack==1/length(TemplateInfo.meansignal)
            % Use mean signal and principal components from TemplateInfo 
            % to (re)calculate layer templates using the appropriate values 
            % of model order and pcPolOrder. 
            if strcmp(Model.type,'PCA')
                Template_all(j) = polyapprox(TemplateInfo.meansignal,...
                    TemplateInfo.pc,Model); % 2014-08-21 14:08
            elseif strcmp(Model.type,'FFTcomp')
                disp('FFTcomp not finalized yet')
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

    %% Processed data for interval:
    % Model array corresponding to species j:
    ModelTemplate = Model;
    ModelTemplate.species = Model.species(j); 
    ModelTemplate.nSpecies = 1;
    ModelTemplate.dstart = Model.manualtemplates(1);
    ModelTemplate.dend = Model.manualtemplates(2);
       
    % Does "Data" contain data for the appropriate interval?
    if Data.depth(1)<=Model.manualtemplates(1) && ...
            Data.depth(end)>=Model.manualtemplates(2)       
        
        % Finalize processing (necessary if floating distances):
        preprocsteps=setpreprocdist(Model.preprocess(:,2),meanLambda);
        [DataProcessed.data, DataProcessed.depth] = makedatafile(Data.data(:,:,j),Data.depth,...
            preprocsteps,Model.derivatives); % No further downsampling or plotting
        
    else
        % Load or, maybe, load and process data for interval.
        % Processed data may not exist from before if the manual layer 
        % template interval is not covered by data section used for layer 
        % counting.

        % Construct processed data file: 
        % All preprocessing steps must be performed. 
        ModelTemplate.preprocess_init = preprocessTemplate(j);
        DataProcessed = constructdatafile(ModelTemplate,manualcounts,Runtype); % 2014-08-23 13:58
    end
    
    %% Calculate layer templates:
    [Template, TemplateInfo] = layerstructure(DataProcessed.data,...
        DataProcessed.depth,manualcounts(:,1),manualcounts(:,3),...
        ModelTemplate,Runtype); % 2014-08-22 15:25
    
    %% Save templates:
    % As data file:
    save([filename '.mat'],'TemplateInfo')
    
    % ...and in array "Template_all":
    Template_all(j)=Template;
    TemplateInfo_all(j)=TemplateInfo;
    
    %% Plot principal components, and save figure:
    if Runtype.plotlevel > 0
        if strcmp(Model.type,'PCA')
            plotprincomp(TemplateInfo,ModelTemplate)
            print([filename '.jpeg'], '-djpeg','-r400')
            % Close plot (maybe):
            if Runtype.plotlevel==1; close(gcf); end
        end
    end
end
end

function plotprincomp(TemplateInfo,Model)
% Plot mean signal and principal components of layers (Model array here 
% only contains one species).
% Mai Winstrup

%% Mean signal:
figure;
subplot(2,1,1)  
x = 1/(2*length(TemplateInfo.meansignal)):1/length(TemplateInfo.meansignal):1;
plot(x,TemplateInfo.meansignal,'-k','linewidth',2)
% Add label and title:
ylabel('Mean signal','fontweight','bold')
if strcmp(Model.icecore,'SyntheticData')
    title(['Species #' Model.species{1} ' (' Model.manCountsName ')'],'fontweight','bold')
else
    title([Model.species{1} ' (' Model.manCountsName ')'],'fontweight','bold')
end

%% Principal components (only plotting the first three components): 
subplot(2,1,2)
plot(x,TemplateInfo.pc(:,1:3),'linewidth',2)
ylabel('PCs','fontweight','bold')
% Add legend:
for i = 1:3; legendname{i,:} = ['PC' num2str(i)]; end
legend(legendname)
legend('boxoff')
end