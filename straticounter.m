function straticounter(sett_icecore)

%% StratiCounter: A layer counting algorithm
% The algorithm is based on the principles of statistical inference of 
% hidden states in semi-Markov processes. States, and their associated 
% confidence intervals, are inferred by the Forward-Backward algorithm and 
% the Viterbi algorithm (optional). The EM (Expectation-Maximization) 
% algorithm is used to find the optimal set of layer parameters for each 
% data batch. Confidence intervals do not account for the uncertainty in 
% estimation of layer parameters.
%
% If (absolute) tiepoints are given, the algorithm is run between these,
% while assuming constant annual layer signals between each pair. If no
% tiepoints, the algorithm is run batch-wise down the core, with a slight 
% overlap between consecutive batches.
%
% See Winstrup (2011) and Winstrup et al. (2012) for further documentation. 
%
% The algorithm was developed for visual stratigraphy data from the NGRIP
% ice core (Winstrup (2011), Winstrup et al. (2012)). It has later been 
% applied to other cores, and extended to parallel analysis of multi-
% parameter data sets (e.g. Vallelonga et al. (2014), Sigl et al. (in prep, 
% 2015)). For testing purposes, it can also be run on synthetic data. 
%
% Developed by Mai Winstrup. 
% Contact: mai@gfy.ku.dk
%
% When using this script, please provide release date of the algorithm, 
% and cite: 
% Winstrup et al., An automated approach for annual layer counting in
% ice cores, Clim. Past. 8, 1881-1895, 2012.
clc; close all;
releasedate = '02-02-2015';

% Copyright (C) 2015  Mai Winstrup
% Files associated with the matchmaker software (matchmaker.m, 
% matchmaker_evaluate.m) is authored and copyrighted by Sune Olander 
% Rasmussen. 
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.

% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
% Public License for more details.

% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

%% Check that settings file exist:
if ~exist(['./Settings/' sett_icecore '.m'],'file')
    disp('Settings file unknown, please correct')
    return
end

%% Add paths to subroutines and settings folders:
addpath(genpath('./Subroutines'))
addpath(genpath('./Settings'))

%% Select how to run the script:
Runtype.develop = 'no';
% In development mode; will run as normal, but output will be put in the
% ./Output/develop folder. Option to run for fewer batches. 
Runtype.reuse = 'yes';
% If yes; use previously processed data and calculated layer templates.
% If no, these are re-calculated. 
Runtype.plotlevel = 1;
% Options: 0: 'none' (no plots), 1: 'info' (few plots), 2: 'debug' (all plots)

% Display info messages if different from standard settings:
if strcmp(Runtype.develop,'yes'); 
    disp('Running in development mode...'); 
end
if strcmp(Runtype.reuse,'no'); 
    disp('Will reprocess data and recalculate layer templates'); 
end
if Runtype.plotlevel>1; 
    disp('Plots will be generated'); 
end

%% Select model settings:
% Import default settings:
Model = defaultsettings; % 2015-01-20

% Use core-specific settings:
run(sett_icecore)

% Add release date:
Model.releasedate = releasedate;

%% Ensure correct format of content in Model:
Model = adjustmodel(Model); % 2015-01-21

%% Make output folders:
[outputdir, outputdir0, runID] = makeoutputfolder(Model, Runtype); % 2015-01-21

%% Load data and manual layer counts:
if strcmp(Model.icecore,'SyntheticData')
    % Check Model.SynthData, and convert mean signal etc. to polynomial 
    % approximations:
    Model.SynthData = checksyntheticmodel(Model,outputdir,Runtype.plotlevel); % 2014-08-23

    % Construct synthetic data: 
    [Data, manualcounts, Model] = makesyntheticdata(Model,Runtype,outputdir); % 2014-08-21
    
else
    % Load manually counted annual layers:
    [manualcounts, meanLambda] = loadlayercounts(Model,[Model.dstart Model.dend]); % 2015-01-21
    
    % Check format, and convert ages to ageUnitOut:
    [manualcounts, Model] = adjustmanualcounts(manualcounts,Model); % 2014-08-23
    
    % Check preprocessing distances (Model.preprocess, Model.dx) relative 
    % to manual layer thicknesses over interval:
    checkpreprocessdist(Model,manualcounts); % 2014-10-01
    
    % Load and preprocess data files:
    [Data, Model] = constructdatafile(Model,manualcounts,Runtype); % 2014-10-01
end

% If no manual counts are known, an empty array should be provided.
% Similarly if information from manual counts should be disregarded:
% manualcounts = [];

%% Save an updated version of "Model" in output folder:
save([outputdir '/Model.mat'],'Model')

%% Initial layer templates and layer parameters:
% These are based on manual layer counts in the data in the depth interval 
% given in "Model.manualtemplates". Layer templates are computed solely 
% based on the data series, not their derivatives.
if isempty(Model.manualtemplates)
    disp('Could allow for using e.g. a sinusoidal layer template. Not implemented.')
else
    [Template0,Template0Info] = ...
        constructmanualtemplates(Data,Model,outputdir,Runtype); % 2014-10-01 13:21
end

% Plot and save figure of initial layer templates:
if Runtype.plotlevel > 0
    % Color of mean signal and PCs:
    color = [0.5 0.5 0.5; 0 0 1; 0 1 0; 1 0 0];
    filename = [outputdir '/layertemplates.jpeg'];        
    hfig_basis = plotlayertemplates(Template0,...
            Template0Info,Model,nan,color,filename); % 2014-08-21 12:56 
    % Close figure?
    if Runtype.plotlevel==1; close(hfig_basis); end
end

% Layer parameters:
if isfield(Model,'Layerpar0')
    % Using a prescribed set of initial parameters:
    Layerpar0 = Model.Layerpar0;
else
    % Layer parameters based on manual layer counts the above calculated 
    % initial layer templates. Parameters are calculated based on data and 
    % their derivatives. 
    Layerpar0 = constructmanualpar(Data,Template0,Model,outputdir,Runtype);
end

% If no parameter updates: Use the values from manual counts?
noupdates = strcmp(Model.update,'none');
index = find(noupdates==1);
names = {'my','sigma','par','cov','nvar'};
for i= 1:length(index)
    disp(['The value of ' names{index(i)} ' will be held constant.'])
    disp('Value corresponding to manual counts is: ');
    disp(num2str(eval(['Layerpar0.' names{index(i)}])))
    reply = input('Use this value? ','s');
    if strcmp(reply,''); disp('yes'); end
    
    if ~ismember(reply,{'yes','y',''})
        value = input('Select new value: ');
        % Check for correct format:
        while ~isequal(size(value),size(eval(['Layerpar0.' names{index(i)}])))
            value = input('Incorrect format, please correct: ');
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

% Save Layerpar0 as output:
save([outputdir '/Layerpar0'],'Layerpar0')

%% Set initial conditions, and initialize arrays: 
[nBatch,batchStart,Layer0,Template,Prior,Layerpar,dDxLambda,logPobs,...
    relweight,Result] = setinitialconditions(Data,Model,manualcounts,...
    meanLambda,Template0,Layerpar0,Runtype); % 2014-10-08 20:41

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BATCHWISE DETECTION OF ANNUAL LAYERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Runtype.develop,'yes')
    nBatch0 = input('Number of batches? (use inf or press return to run for all data) ');
    if isempty(nBatch0); nBatch0=inf; end
    nBatch = min(nBatch,nBatch0);
end

disp('Algorithm is running, please be patient...')
logPobs_alldata = 0; 
iBatch = 0;
   
while iBatch < nBatch
    %% Batch number:
    iBatch = iBatch+1;
    if mod(iBatch,1)==0
        disp(['Batch ' num2str(iBatch) ': ' num2str(Data.depth(batchStart(iBatch))) 'm'])
    end
    
    %% 1: Select data corresponding to current batch
    % And give an upperbound estimate of the number of layers in batch,
    % based on the batch size.
    if isempty(Model.tiepoints)
        % Select data section containing approximately "nLayerBatch":
        meanLambda = exp(Prior(iBatch).m+Prior(iBatch).sigma^2/2); %[m]
        % Length of data sequence (in pixels): 
        batchLength = round(Model.nLayerBatch*meanLambda/Model.dx); 
        batchEnd = min(batchStart(iBatch)+batchLength-1,length(Data.depth));            
        % Upper-bound estimate of years (incl. uncertainties) in batch:
        nLayerMax = round(2*Model.nLayerBatch);
        
        % In case there is less left than a full batch of data, the current 
        % batch is extended to include these extra pixels:
        if batchEnd+batchLength > length(Data.depth)
            % Extended batch size:
            batchEnd = length(Data.depth);
            batchLength = batchEnd-batchStart(iBatch)+1;
            
            % Maximum number of layers in batch is then: 
            nLayer = (Data.depth(end)-Data.depth(batchStart(iBatch)))/meanLambda;
            nLayerMax = round(2*nLayer);

            % Do not use any overlap section for the last data batch:
            Model.batchOverlap = 0;
            disp(['Stopping at batch ' num2str(iBatch) ', ' num2str(Data.depth(batchEnd)) 'm'])
        end
        
    else
        % Using data sequence up to (and including) the next tiepoint:
        batchEnd = interp1(Data.depth,1:length(Data.depth),Model.tiepoints(iBatch+1,1),'nearest');
        batchLength = batchEnd-batchStart(iBatch)+1; 
        % Number of layers in interval, including start/end layers
        % (the actual number is also the maximum number of layers).
        nLayerMax = abs(Model.tiepoints(iBatch+1,2)-Model.tiepoints(iBatch,2))+1;
        % Mean layer thickness:
        meanLambda = exp(Prior(iBatch).m+Prior(iBatch).sigma^2/2); %[m]
    end
    
    % Display warning if the data section may be too small for reliable
    % layer parameter re-estimation (may e.g. occur when using closely 
    % spaced tiepoints):
    if nLayerMax <= 50
        disp(['Note: nLayerMax for batch is small (' num2str(nLayerMax) ')'])
    end
    
    % Depth scale for batch:
    depth_batch = Data.depth(batchStart(iBatch):batchEnd);
              
    %% Process batch data (if using floating preprocessing distance): 
    % Floating distances are based on estimated mean layer thickness 
    % for batch: 

    %   [depth_batch,data_batch] = finalizeprocessing(Data.depth,Data.data,meanLambda,interval,Model,loglevel);
%         preprocess_batch = setfloatingdist(Model.preprocess_rep,meanLambda); % 2014-07-16 16:18
%         
%         data_batch = nan(batchLength,3,Model.nSpecies);
%         for j = 1:Model.nSpecies 
%             % Using a slightly extended area: 
%             Lext = max(cell2mat(preprocess_batch{j}(:,2)))/Model.dx;
%             if isempty(Lext); Lext = 0; end
%             px_start_ext = max(1,px_start(iBatch)-Lext);
%             batchend_ext = min(batchend+Lext,length(Data.depth));
%             
%             rawdata_batch_ext = Data.data(px_start_ext:batchend_ext,1,j);
%             depth_batch_ext = Data.depth(px_start_ext:batchend_ext);
%             
%             procdata_batch_ext = processdata(rawdata_batch_ext,depth_batch_ext,...
%                 Model.species{j},preprocess_batch{j},Model); %2014-07-16 16:18
%             procdata_batch = procdata_batch_ext(1+px_start(iBatch)-px_start_ext:length(depth_batch_ext)-batchend_ext+batchend); 
%             
%             [slope,dslope] = calculateslope(procdata_batch,Model.slopeorder,Model.slopedist,0); %2014-07-16 14:05
%             data_batch(:,:,j) = [procdata_batch, slope, dslope]; % men vi vil jo ikke bruge dem alle!!!
%         end

    data_batch = Data.data(batchStart(iBatch):batchEnd,:,:);
        
    %% 2+3: Iterations over layer template and layer parameters 
    % The layer parameters are first optimized according to the initial 
    % template. Then the layer template is optimized, and layer parameters 
    % are re-optimized for this new template, etc. 

    % Initial layer template for batch:
    % The same initial template is used for all batches:
    Template(:,iBatch,1)=Template(:,1,1);

    % Iterations over layer template:
    for iTemplateBatch = 1:Model.nTemplateBatch
        
        % Iterations over layer parameters:
        iIter = 1;
        while iIter <= Model.nIter
            if nBatch < 5; disp(['Iteration #' num2str(iIter)]); end

            %% 2a: Select model parameters
            if iIter == 1
                % Initial estimate of the selected model parameters is 
                % based on the maximum likelihood of the prior:
                Layerpar(iBatch,iTemplateBatch,1).my = Prior(iBatch).m;
                Layerpar(iBatch,iTemplateBatch,1).sigma = Prior(iBatch).sigma;
                Layerpar(iBatch,iTemplateBatch,1).par = Prior(iBatch).u;
                Layerpar(iBatch,iTemplateBatch,1).cov = Prior(iBatch).cov;
                Layerpar(iBatch,iTemplateBatch,1).nvar = Prior(iBatch).nvar; 
            else 
                % Update model parameters:
                Layerpar(iBatch,iTemplateBatch,iIter) = Layerpar_new;
            end
            
            %% 2b: Run layer detection algorithm
            [Layerpos_new, FBprob_new, ExpVal_new, logPobs_new, d, pd, logb_new, bweight] = ...
                layerdetection(data_batch,Template(:,iBatch,iTemplateBatch),...
                Layerpar(iBatch,iTemplateBatch,iIter),Layer0(iBatch),...
                nLayerMax,batchLength,Model,Runtype); % 2014-10-10 16:37
            
            %% 2c: New estimates for layer parameters
            [Layerpar_new, relweight(iBatch,iTemplateBatch,iIter)] = ...
                updatelayerpar(ExpVal_new,FBprob_new,Prior(iBatch),...
                Layerpar(iBatch,iTemplateBatch,iIter),d,pd,batchLength,...
                logb_new,bweight,Model);
            Layerpar(iBatch,iTemplateBatch,iIter) = Layerpar_new;
            
            % In case of Bayesian estimates, logPobs is re-calculated
            % to account for prior knowledge of layer parameters:
            logPobs_new = calclogPobs(logPobs_new,Layerpar_new,Prior,Model);
            logPobs(iBatch,iTemplateBatch,iIter+1,:) = logPobs_new;

             % Save resulting set of layer positions, potentially both 
             % from Forward-Backward and Viterbi algorithms: 
    %        Layerpos(iBatch,iTemplateBatch,iIter) = Layerpos_new;
            
            %% 2d: Stopping criteria reached?
            % Using logPobs from FB algorithm to test for convergence:
            if abs(logPobs_new(1)-logPobs(iBatch,iTemplateBatch,iIter,1))>Model.eps && ...
                    iIter<=Model.nIter;
                iIter = iIter+1;          
            else
                % Actual number of iterations for this batch and template:
                nIter = iIter;
                if abs(logPobs_new(1)-logPobs(iBatch,iTemplateBatch,iIter,1))<Model.eps
                    flag = 1; break % Convergence reached
                else flag = 2; break % Maximum number of iterations reached
                end
            end
       end
           
       %% 2e: Check that log(Pobs) is always growing (as it should)
       if Runtype.plotlevel>1
           if ~exist('hfig_logPobs','var'); hfig_logPobs = figure; end
           if iTemplateBatch == 1; clf(hfig_logPobs); end
           figure(hfig_logPobs)
           plot(squeeze(logPobs(iBatch,iTemplateBatch,2:end,1)),'.-k'); 
           hold on
           title(['Batch: ' num2str(iBatch)],'fontweight','bold')
           for i = 1:iTemplateBatch; legendname{i} = ['template iteration: ' num2str(i)]; end
           legend(legendname,'location','best')
           ylabel('log(P_{obs})')
       end
       
       %% 3a: Find new layertemplates
       layercounts_new = Data.depth(Layerpos_new.fb+batchStart(iBatch)-1.5)+Model.dx/2;  
       Model_template = Model;
       Model_template.layerCharInterval = [depth_batch(1) depth_batch(end)];
       [Template_new, TemplateInfo_new] = layerstructure(data_batch,...
           depth_batch,layercounts_new,[],Model_template,Runtype); % 2014-08-21 23:22
       
       % Plot new layer templates (and compare to original ones):
       %color = [1 0 1; 1 1 0; 1 0.5 0.5; 0 0 1];
       %filename = [outputdir '/layertemplates_new'];
       %plotlayertemplates(Template_new,TemplateInfo_new,Model,hfig_basis,color,filename);       
       %if Runtype.plotlevel < 1; close(hfig_basis2); end
       
       %% 3b: Update layer templates: 
       Template(:,iBatch,iTemplateBatch+1)=Template_new; 
       % (But first template for each batch is always the same)

%        % With updates:
%        Template(:,iBatch,iTemplateBatch+1)=Template_new;
%        % Only updating the mean trajectory: 
%        for j = 1:Model.nSpecies
%            Template(j,iBatch,iTemplateBatch+1) = Template(j,1,1);
%            Template(j,iBatch,iTemplateBatch+1).mean = Template_new(j).mean;
%        end
    end

    %% 4: When all iterations have converged: 
    % Results for batch:
    [Result(iBatch),Layer0(iBatch+1),batchStart(iBatch+1)] = ...
        resultsforbatch(depth_batch,FBprob_new,Layerpos_new,...
        Layer0(iBatch),d,pd,logb_new,meanLambda,batchStart(iBatch),dDxLambda,...
        iIter,Model,Runtype);
    
    % Calculate log(P_obs) for entire data series, up and including this one:
    logPobs_alldata = logPobs_alldata + logPobs(iBatch,iTemplateBatch,iIter,1); 
    % (up to tau in this batch) - no to T? if so, there is an overlap...
    
    % Prior for next batch:
    Prior(iBatch+1) = updatepriors(Prior(iBatch),Layerpar_new,Model);

    %% 5: Visualize results and compare to manual layer counts:
    if mod(iBatch,5)==0
        % Combine results from batches, and convert from layers to ages:
        [Layerpos_prelim,~,~,timescale_prelim,timescale1yr_prelim,~,...
             markerConf_prelim,lambda_prelim] = ...
             combinebatches(Result(1:iBatch),manualcounts,Model); % 2014-08-12 21:58
            
        % Save preliminary output:
        save([outputdir '/timescale1yr_prelim'],'timescale1yr_prelim')
        save([outputdir '/markerhorizons_prelim'],'markerConf_prelim')
           
        % Plot preliminary timescale and mean layer thicknesses:
        if Runtype.plotlevel>0
            % Timescale:
            if exist('hfig_timescale','var'); close(hfig_timescale); end
            filename = [outputdir '/timescale_prelim.jpeg'];
            hfig_timescale = showtimescale(timescale_prelim,...
                timescale1yr_prelim,Layerpos_prelim,manualcounts,...
                Data.depth(batchStart(1:iBatch)),Model,filename); % 2014-08-21 23:32
            
            % Mean layer thicknesses:
            if exist('hfig_lambda','var');
                for i = 1:length(Model.dxLambda)
                    if ~isempty(hfig_lambda(i)); close(hfig_lambda(i)); end
                end
            end
            hfig_lambda = gobjects(length(Model.dxLambda),1); % Initialize handles
            for idx = 1:length(Model.dxLambda)                
                filename = [outputdir '/lambda_' num2str(Model.dxLambda(idx)) 'm_prelim.jpeg'];
                hfig_lambda(idx) = showlambda(lambda_prelim{idx},...
                   Layerpos_prelim,timescale1yr_prelim,manualcounts,Model,filename); % 2014-08-15 11:12               
            end
        end
    end
 
    %% Reached end of data record?
    if isequal(batchEnd,length(Data.depth)); break; end
    % Or, alternatively, the last tiepoint:
    if iBatch==size(Model.tiepoints,1)-1; break; end

    % If not: Expand initialized matrices?
    if iBatch==nBatch
        % Remaining data interval:
        depthinterval = Model.dend-Data.depth(batchEnd);
        
        % Estimate remaining number of batches:
        meanLambda = exp(Prior(iBatch).m+Prior(iBatch).sigma^2/2); %[m]
        nBatchRest = ceil(1.1*depthinterval/(meanLambda*Model.nLayerBatch));
        % New value of nBatch:
        nBatchNew = nBatch+nBatchRest;
        
        % If running in develop mode, we may have provided a maximum batch 
        % number:
        if strcmp(Runtype.develop,'yes')
            nBatchNew = min(nBatchNew,nBatch0);
            nBatchRest = nBatchNew-nBatch;
        end
        nBatch = nBatchNew;
        
        % Expand matrices:
        if nBatchRest > 0
            [batchStartRest,Layer0Rest,TemplateRest,PriorRest,...
                LayerparRest,logPobsRest,relweightRest,ResultRest]=...
                initializematrices(nBatchRest,Model);
            batchStart = [batchStart; batchStartRest];
            Layer0 = [Layer0, Layer0Rest];
            Template = [Template, TemplateRest];
            Prior = [Prior; PriorRest];
            Layerpar = [Layerpar; LayerparRest];
            logPobs = [logPobs; logPobsRest];
            relweight = [relweight; relweightRest];
            Result = [Result, ResultRest];
        else
            break
        end
    end
end 

%% Combine results from batches: 
% Actual number of batches:
nBatch = iBatch;
% Remove unused parts of initialized matrices:
batchStart = batchStart(1:nBatch);
Result = Result(1:nBatch);
 
% Combine batches:
[Layerpos,LayerDist,centralEst,timescale,timescale1yr,markerProb,...
    markerConf,lambdaResults] = combinebatches(Result,manualcounts,Model); % 2014-08-12 21:58 

%% Save output:
% Timescale:
save([outputdir '/timescale.mat'],...
    'timescale','timescale1yr','Layerpos','LayerDist','centralEst','Model')
% Save timescale1yr as textfile with metadata:
filename = [outputdir '/' Model.icecore '_timescale1yr.txt'];
savetimescaleastxt(timescale1yr,filename,Model) % 2014-08-14 14:41

% Confidence interval for marker horizons:
if ~isempty(Model.dMarker)
    save([outputdir '/markerhorizons.mat'],...
        'markerProb','markerConf','Model')
    % As txt file:
    filename = [outputdir '/' Model.icecore '_markerconf'];
    savemarkerastxt(markerConf,filename,Model) % 2014-10-22 18:03
end

% Layer thicknesses:
if ~isempty(Model.dxLambda)
    save([outputdir '/lambda.mat'],'lambdaResults','Model')
end
   
% Save run number:
save([outputdir0 '/runID.mat'],'runID')
 
%% Plot timescale and mean layer thicknesses:
if Runtype.plotlevel>0 
    % Timescale:
    if exist('hfig_timescale','var'); close(hfig_timescale); end
    filename = [outputdir '/timescale.jpeg'];
    hfig_timescale = showtimescale(timescale,timescale1yr,Layerpos,...
        manualcounts,Data.depth(batchStart),Model,filename); % 2014-08-21 23:32
    
    % Mean layer thicknesses:
    if exist('hfig_lambda','var');
        for i = 1:length(Model.dxLambda)
            if ~isempty(hfig_lambda(i)); close(hfig_lambda(i)); end
        end
    end
    hfig_lambda = gobjects(length(Model.dxLambda),1); % Initialize handles
    for idx = 1:length(Model.dxLambda)
        filename = [outputdir '/lambda_' num2str(Model.dxLambda(idx)) 'm.jpeg'];
        hfig_lambda(idx) = showlambda(lambdaResults{idx},Layerpos,...
            timescale1yr,manualcounts,Model,filename); % 2014-08-15 11:12
    end
end
   
%% Calculate new layer templates:
% Based on data in complete interval.
[TemplateNew,TemplateInfoNew]=layerstructure(data_batch,depth_batch,...
    Layerpos.combined,[],Model,Runtype); % 2014-08-21 23:22

%% Clean-up: Remove preliminary datafiles and figures
% Data files:
filename = [outputdir '/timescale1yr_prelim.mat'];
if exist(filename,'file'); delete(filename); end
filename = [outputdir '/markerhorizons_prelim.mat'];
if exist(filename,'file'); delete(filename); end
% Figures:
filename = [outputdir '/timescale_prelim.jpeg'];
if exist(filename,'file'); delete(filename); end
for idx = 1:length(Model.dxLambda)
    filename = [outputdir '/lambda_' num2str(Model.dxLambda(idx)) 'm_prelim.jpeg'];
    if exist(filename,'file'); delete(filename); end
end

%% Display name of output directory:
disp(['Output directory: ' outputdir])

%% Show results in matchmaker:
checkinmatchmaker(outputdir,Model);