function straticounter(varargin)

%% STRATICOUNTER: A layer counting algorithm
% STRATICOUNTER(settings_file_name) loads settings from the named file in
%   "settings_file_name". 
%
% STRATICOUNTER(settings_path, output_path) loads settings from a ".mat"
%   file whose path and name is given in "settings_path". The output of the
%   function will be saved in the "output_path".
%
% INPUTS:
%   settings_file_name (string): Name of the settings file (.m) to load.
%       It is assumed that the file is located in a local folder
%       called "Settings".
%   settings_path (string): Path to the file (.mat) from which the
%       settings should be loaded.
%   output_path (string): Path to the folder where the output of the
%       function should be written.
%
% OUTPUTS:
%   When providing only 1 argument, the output files will be placed in a
%   local folder named "Output"
%
%   When providing 2 arguments, the output files will be saved in the
%   folder specified by the "output_path" argument.
%
% DETAILS:
% The algorithm is based on the principles of statistical inference of
% hidden states in semi-Markov processes. States, and their associated
% confidence intervals, are inferred by the Forward-Backward algorithm.
% The EM (Expectation-Maximization) algorithm is used to find the optimal
% set of layer parameters for each data batch. Confidence intervals do not
% account for the uncertainty in estimation of layer parameters.
%
% If (absolute) tiepoints are given, the algorithm is run between these,
% while assuming constant annual layer signals between each pair. If no
% tiepoints, the algorithm is run batch-wise down the core, with a slight
% overlap between consecutive batches.
%
% The algorithm was developed for visual stratigraphy data from the NGRIP
% ice core (Winstrup (2011), Winstrup et al. (2012)). It has later been
% applied to other cores, and extended to parallel analysis of multi-
% parameter data sets (e.g. Vallelonga et al. (2014), Sigl et al. (in prep,
% 2015)). For testing purposes, it can also be run on synthetic data.
%
% See Winstrup (2011) and Winstrup et al. (2012) for further documentation.
%
% When using this script, please provide release date of the algorithm,
% and cite:
% Winstrup et al., An automated approach for annual layer counting in
% ice cores, Clim. Past. 8, 1881-1895, 2012.
%
%% Copyright (C) 2015  Mai Winstrup
% Files associated with the matchmaker software (matchmaker.m,
% matchmaker_evaluate.m) is authored and copyrighted by Sune Olander
% Rasmussen.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

%% Release date:
releasedate = '07-07-2015';

%% Paths to subroutine and settings folders:
if ~isdeployed % shows whether this is a compiled instance of the code
    % Header for current run:
    disp('==========================================')
    disp(varargin{1})
    disp('==========================================')
    % close all
    addpath(genpath('./Subroutines'))
    addpath(genpath('./Settings'))
end

%% Select model settings:
% Import default settings:
Model = defaultsettings();
% Add release date:
Model.releasedate = releasedate;

%% Core-specific settings:
% Make generic message to be used when incorrect number of inputs are used:
vararg_err = 'This function accepts a maximum of two (2) input arguments';

% Use 'nargin' to see the number of input arguments, and select appropriate 
% behavior:
if nargin == 1    
    % Check that settings file exists:
    if ~exist(['Settings/' varargin{1}],'file')
        error('Settings file unknown, please correct')
    end
    % Import settings:
    run(varargin{1});
elseif nargin ==2
    run(varargin{1});
else
    error(vararg_err);
end

%% Select how to run the script:
if nargin == 1
    Runtype.develop = 'no';
    % In development mode; will run as normal, but output will be put in 
    % the ./Output/develop folder. Option to run for only a few batches.
    Runtype.reuse = 'yes';
    % If yes; use previously processed data and calculated layer templates.
    % If no, these are re-calculated.
    Runtype.outdir = '';
    % Provide path if wanting to use a custom output directory.
    
elseif nargin == 2
    Runtype.develop = 'no';
    Runtype.reuse = 'no';
    Runtype.outdir = varargin{2};
else
    error(vararg_err)
end

% Plotting:
% Options: 0: 'none' (no plots), 1: 'info' (few plots), 2: 'debug' (all plots)
if isdeployed
    % Force to no plots when run as compiled library
    Runtype.plotlevel = 0;
else
    Runtype.plotlevel = 1;
end

% Display info messages if different from standard settings:
if strcmp(Runtype.develop,'yes');
    disp('Running in development mode...');
end
if strcmp(Runtype.reuse,'no');
    disp('Will reprocess data and recalculate layer templates');
end
if Runtype.plotlevel>1;
    disp('Many plots will be generated');
end

%% Ensure correct format of content in Model:
Model = adjustmodel(Model);

%% Make output folders:
[outputdir, outputdir0, runID] = makeoutputfolder(Model, Runtype);

%% Load data and manual layer counts:
if strcmp(Model.icecore,'SyntheticData')
    % Construct synthetic data:
    [Data, manualcounts, Model] = makesyntheticdata(Model,Runtype,outputdir);

else
    % Load manually-counted annual layers for interval:
    [manualcounts, meanLambda] = loadlayercounts(Model,[Model.dstart Model.dend]);
    % If no manual counts are known, an empty array is provided, and
    % meanLambda is estimated by user:
    while isnan(meanLambda) || meanLambda <= 0
        disp(['Rough estimate for mean layer thickness (in m) '...
            'over layer counting interval'])
        meanLambda = input(['(used for finding data gabs and '...
            'estimating value of nBatch): ']);
    end

    % Check format, and convert ages to ageUnitOut:
    [manualcounts, Model] = adjustmanualcounts(manualcounts,Model);

    % Check preprocessing distances relative to manual layer thicknesses
    % over interval:
    checkpreprocdist(Model,manualcounts);

    % Load and preprocess data files:
    [Data, Model] = loadormakedatafile(Model,manualcounts,Runtype);

    % Check for long sections without data:
    % Long sections are in this context corresponding to 20 mean layer
    % thicknesses without much data:
    sectionswithoutdata(Data,20*meanLambda,Model.species);
end

%% Initial layer templates and layer parameters:
% These are based on manual layer counts in the data in the depth interval
% given in "Model.manualtemplates". Layer templates are computed solely
% based on the data series, not their derivatives.
if isempty(Model.manualtemplates)
    disp('Could allow for using e.g. a sinusoidal layer template. Not implemented.')
else
    [Template0,Template0Info,Model] = ...
        constructmanualtemplates(Data,Model,outputdir,Runtype);
end

% Plot and save figure of initial layer templates:
if Runtype.plotlevel > 0
    % Color of mean signal and PCs:
    color = [0.5 0.5 0.5; 0 0 1; 0 1 0; 1 0 0];
    filename = [outputdir '/layertemplates.jpeg'];
    hfig_template = plotlayertemplates(Template0,...
            Template0Info,Model,nan,color(1:Model.order+1,:),filename);
end

% Layer parameters:
if isfield(Model,'Layerpar0')
    % Using a prescribed set of initial parameters:
    Layerpar0 = Model.Layerpar0;
else
    % Layer parameters based on manual layer counts and the layer templates
    % calculated above. Parameters are calculated based on data and their
    % derivatives.
    Layerpar0 = constructmanualpar(Data,Template0,Model,outputdir,Runtype);
end

%% Save an updated version of "Model" in output folder:
save([outputdir '/Model.mat'],'Model')

%% Set initial conditions, and initialize arrays:
[nBatch,batchStart,Layer0,Template,Prior,Layerpar,dDxLambda,logPobs,...
    logPobsNorm,relweight,Result] = setinitialconditions(Data,Model,...
    manualcounts,meanLambda,Template0,Layerpar0);

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
data_final = [];

while iBatch < nBatch    
    %% Batch number:
    iBatch = iBatch+1;
    if nBatch<10; dispBatch = 1; else dispBatch = 5; end
    if mod(iBatch,dispBatch)==0
        disp(['Batch ' num2str(iBatch) ': ' ...
            num2str(Data.depth(batchStart(iBatch))) 'm'])
    end

    %% 1a: Select data corresponding to current batch:
    % And provide an upper-bound estimate of the number of layers in batch.
    if isempty(Model.tiepoints)
        % Situation A: No tiepoints.
        % We select a section of appropriate length. This section may be
        % larger than nLayerBatch if close to the end of data series.

        % Select data section containing approximately "nLayerBatch":
        meanLambda = exp(Prior(iBatch).m+Prior(iBatch).sigma^2/2); %[m]
        % Length of data sequence (in pixels):
        batchLength = round(Model.nLayerBatch*meanLambda/Model.dx);
        batchEnd = min(batchStart(iBatch)+batchLength-1,length(Data.depth));
        % Upper-bound estimate of years (incl. uncertainties) in batch:
        nLayerMax = round(1.3*Model.nLayerBatch);
        % If indications that this value should be larger, it will be
        % increased with each batch iteration.

        % In case there is less left than a full batch of data, the current
        % batch is extended to include these extra pixels:
        if batchEnd+batchLength > length(Data.depth)
            % Extended batch size:
            batchEnd = length(Data.depth);
            batchLength = batchEnd-batchStart(iBatch)+1;

            % Maximum number of layers in batch is then:
            nLayerMean = (Data.depth(end)-Data.depth(batchStart(iBatch)))/meanLambda;
            nLayerMax = round(1.3*nLayerMean);            
            
            % Do not use any overlap section for the last data batch:
            Model.batchOverlap = 0;
            disp(['Stopping at batch ' num2str(iBatch) ', ' ...
                num2str(Data.depth(batchEnd)) 'm'])
        end

    else
        % Situation B: Tiepoints. Using sections between tiepoints.
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
        disp(['Note: Batch contains a limited number of layers (nLayerMax = ' ...
            num2str(nLayerMax) ')'])
        disp(['For reliable layer parameter re-estimation it is recommended '...
            'to use a larger batch size.'])
    end

    % Depth scale for batch:
    depth_batch = Data.depth(batchStart(iBatch):batchEnd);

    %% 1b: Preprocess batch data:
    % Data in batch is preprocessed relative to the mean layer thickness as
    % estimated from previous batch.
    % Preprocessing specs:
    [preprocsteps,preprocdist] = ...
        setpreprocdist(Model.preprocsteps(:,2),meanLambda);

    % Using an extended section around batch (extension depends on
    % maximum processing distance):
    L = max(cell2mat(preprocdist));
    if isempty(L); L = 0; end
    istart = find(Data.depth>=depth_batch(1)-L/2,1,'first');
    iend = find(Data.depth<=depth_batch(end)+L/2,1,'last');
    data_in = Data.data(istart:iend,:,:);
    depth_in = Data.depth(istart:iend);

    % Preprocess batch data (keep depth scale):
    data_out = makedatafile(data_in,depth_in,preprocsteps,Model.derivatives,depth_in);
    % Save resulting data record:
    data_final = [data_final; data_out];

    % Remove extended part of data:
    data_batch = data_out(batchStart(iBatch)-istart+1:batchEnd-istart+1,:,:);

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
            if nBatch < 5; disp(['  > Iteration #' num2str(iIter)]); end

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

            %% 2b: Run layer detection algorithm:
            [Layerpos_new, FBprob_new, ExpVal_new, logPobs_new, d, pd, ...
                logb_new, nLayerMax_new] = ...
                layerdetection(data_batch,Template(:,iBatch,iTemplateBatch),...
                Layerpar(iBatch,iTemplateBatch,iIter),Layer0(iBatch),...
                nLayerMax,Model,Runtype.plotlevel);
            
            %% 2c: New estimates for layer parameters:
            [Layerpar_new, relweight_new] = ...
                updatelayerpar(ExpVal_new,FBprob_new,Prior(iBatch),...
                Layerpar(iBatch,iTemplateBatch,iIter),d,pd,logb_new,...
                Model);

            % log(P_obs): In case of Bayesian estimates, this value is
            % re-calculated to account for prior probability of obtained
            % layer parameters:
            logPobs_new = calclogPobs(logPobs_new,Layerpar_new,Prior,Model);

            %% 2d: Save results for iteration:
            % New set of layer parameters:
            Layerpar(iBatch,iTemplateBatch,iIter+1) = Layerpar_new;
            % Relative weighting factor:
            relweight(iBatch,iTemplateBatch,iIter+1) = relweight_new;
            % Value of log(P_obs):
            logPobs(iBatch,iTemplateBatch,iIter+1) = logPobs_new;

            % New estimate for maximum numbers of layer in batch:
            nLayerMax = nLayerMax_new;                
            
            %% 2e: Stopping criteria reached?
            % Using log(P_obs) to test for convergence:
            if logPobs(iBatch,iTemplateBatch,iIter+1)-...
                    logPobs(iBatch,iTemplateBatch,iIter)>Model.eps && ...
                    iIter<Model.nIter;
                iIter = iIter+1;
                
            else
                % Convergence or maximum iteration value has been reached.
                % Number of iterations performed:
                nIter = iIter;
                if abs(logPobs(iBatch,iTemplateBatch,nIter)-...
                        logPobs(iBatch,iTemplateBatch,nIter+1))<=Model.eps
                    flag = 1; break % Convergence reached
                else flag = 2; break % Maximum number of iterations reached
                end
            end

       end

       %% 2f: Check that log(Pobs) is always growing (as it should)
       if Runtype.plotlevel>1 && Model.nIter > 1
           if iTemplateBatch==1;
              if iBatch == 1;
                  hfig_logPobs = figure;
              else clf(hfig_logPobs); % Clear figure from previous batch
              end
           end
           figure(hfig_logPobs)
           plot(squeeze(logPobs(iBatch,iTemplateBatch,2:end)),'.-',...
               'color',[1 1 1]*(iTemplateBatch-1)/Model.nTemplateBatch);
           hold on
           title(['Batch: ' num2str(iBatch)],'fontweight','bold')
           xlabel('Iteration number')
           ylabel('log(P_{obs})')
           if iTemplateBatch > 1
               for i = 1:iTemplateBatch;
                   legendname{i} = ['Template iteration: ' num2str(i)];
               end
               legend(legendname,'location','best')
           end
       end

       %% 3: Improve the layertemplates:
       layercounts_new = Data.depth(Layerpos_new.fb+batchStart(iBatch)-1.5)+Model.dx/2;
       ModelTemplate = Model;
       ModelTemplate.layerCharInterval = [depth_batch(1) depth_batch(end)];
       [Template_new, TemplateInfo_new] = layerstructure(data_batch,...
           depth_batch,layercounts_new,[],ModelTemplate,Runtype);
       % Save the new layer templates:
       Template(:,iBatch,iTemplateBatch+1)=Template_new;

       % Plot new layer templates (and compare to original ones):
       if Runtype.plotlevel > 1
           color = [1 0 1; 1 1 0; 1 0.5 0.5; 0 0 1];
           filename = [outputdir '/layertemplates_new'];
           plotlayertemplates(Template_new,TemplateInfo_new,Model,...
               hfig_template,color,filename);
       end

%        % Possibility for updating the layer templates:
%        % (Currently, the initial template for each batch is always the same)
%        Template(:,iBatch,iTemplateBatch+1)=Template_new;
%        % Or: Only updating the mean trajectory:
%        for j = 1:Model.nSpecies
%            Template(j,iBatch,iTemplateBatch+1) = Template(j,1,1);
%            Template(j,iBatch,iTemplateBatch+1).mean = Template_new(j).mean;
%        end
    end

    %% 4: When all iterations have converged:
    % Results for batch: layer probability distribution, the most likely
    % layer boundaries, and initial conditions for next batch.
    [Result(iBatch),Layer0(iBatch+1),batchStart(iBatch+1)] = ...
        resultsforbatch(depth_batch,FBprob_new,Layerpos_new,...
        Layer0(iBatch),d,pd,logb_new,meanLambda,batchStart(iBatch),...
        dDxLambda,iIter,Model,Runtype.plotlevel);

    % Calculate log(P_obs) for entire data series, up and including this one:
    logPobs_alldata = logPobs_alldata + logPobs(iBatch,iTemplateBatch,iIter);
    % (OBS: Does this go up to tau in this batch, or to T? If so, the
    % probabilities for ending sections are counted twice.)

    % logPobs normalized to the number of data point in batch (the full
    % batch, since logPobs includes probability section after tau):
    logPobsNorm(iBatch) = logPobs(iBatch,iTemplateBatch,iIter)/batchLength;

    % Prior for next batch:
    Prior(iBatch+1) = updatepriors(Prior(iBatch),Layerpar_new,Model);

    %% 5: Visualize results and compare to manual layer counts:
    if mod(iBatch,5)==0
        % Combine results from batches, and convert from layers to ages:
        [~,~,~,timescale_prelim,timescale1yr_prelim,~,...
             markerConf_prelim,lambda_prelim] = ...
             combinebatches(Result(1:iBatch),manualcounts,Model);
         
        % Save preliminary output:
        save([outputdir '/timescale1yr_prelim'],'timescale1yr_prelim')
        save([outputdir '/markerhorizons_prelim'],'markerConf_prelim')

        % Plot preliminary timescale and mean layer thicknesses:
        if Runtype.plotlevel>0
            % Timescale:
            if exist('hfig_timescale','var'); close(hfig_timescale); end
            filename = [outputdir '/timescale_prelim.jpeg'];
            hfig_timescale = showtimescale(timescale_prelim,...
                timescale1yr_prelim,manualcounts,...
                Data.depth(batchStart(1:iBatch)),Model,filename);

            % Mean layer thicknesses:
            if exist('hfig_lambda','var')
                close(hfig_lambda(isgraphics(hfig_lambda,'figure')));
            end
            hfig_lambda = gobjects(length(Model.dxLambda),1); % Initialize handles
            for idx = 1:length(Model.dxLambda)
                filename = [outputdir '/lambda_' num2str(Model.dxLambda(idx)) 'm_prelim.jpeg'];
                hfig_lambda(idx) = showlambda(lambda_prelim{idx},...
                   timescale1yr_prelim,manualcounts,Model,filename);
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
        nBatchRest = ceil(1.3*depthinterval/(meanLambda*Model.nLayerBatch));
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
            [batchStartRest,Layer0Rest,TemplateRest,PriorRest,LayerparRest,...
                logPobsRest,logPobsNormRest,relweightRest,ResultRest]=...
                initializematrices(nBatchRest,Model);
            batchStart = [batchStart; batchStartRest];
            Layer0 = [Layer0, Layer0Rest];
            Template = [Template, TemplateRest];
            Prior = [Prior; PriorRest];
            Layerpar = [Layerpar; LayerparRest];
            logPobs = [logPobs; logPobsRest];
            logPobsNorm = [logPobsNorm; logPobsNormRest];
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
Result = Result(1:nBatch);
Layerpar = Layerpar(1:nBatch,:,:);
relweight = relweight(1:nBatch,:,:);
logPobs = logPobs(1:nBatch,:,:);
logPobsNorm = logPobsNorm(1:nBatch);
batchStart = batchStart(1:nBatch+1);
% Convert batchStart to depths:
batchStartDepth = Data.depth(batchStart);

% Combine batches:
[Layerpos,LayerProbDist,centralEst,timescale,timescale1yr,markerProb,...
    markerConf,lambdaResults] = combinebatches(Result,manualcounts,Model);

%% Save output:
% Timescale:
save([outputdir '/timescale.mat'],...
    'timescale','timescale1yr','Layerpos','LayerProbDist','centralEst','Model')
% Save timescale1yr as textfile with metadata:
filename = [outputdir '/' Model.icecore '_timescale_' num2str(Model.dstart) ...
    '-' num2str(Model.dend) 'm.txt'];
savetimescaleastxt(timescale1yr,filename,Model)

% Confidence interval for marker horizons:
if ~isempty(Model.dMarker)
    save([outputdir '/markerhorizons.mat'],...
        'markerProb','markerConf','Model')
    % As txt file:
    filename = [outputdir '/' Model.icecore '_markerconf'];
    savemarkerastxt(markerConf,filename,Model)
end

% Layer thicknesses:
if ~isempty(Model.dxLambda)
    save([outputdir '/lambda.mat'],'lambdaResults','Model')
end

% Save all results from iterations (!):
save([outputdir '/results.mat'],'Result','relweight','logPobs','logPobsNorm',...
    'Layerpar','Prior','Template','batchStartDepth','data_final')

% Save run number:
save([outputdir0 '/runID.mat'],'runID')

%% Plot timescale and mean layer thicknesses:
if Runtype.plotlevel>0
    % Timescale:
    if exist('hfig_timescale','var'); close(hfig_timescale); end
    filename = [outputdir '/timescale.jpeg'];
    hfig_timescale = showtimescale(timescale,timescale1yr,...
        manualcounts,Data.depth(batchStart),Model,filename);

    % Mean layer thicknesses:
    if exist('hfig_lambda','var');
        for i = 1:length(Model.dxLambda)
            if ~isempty(hfig_lambda(i)); close(hfig_lambda(i)); end
        end
    end
    hfig_lambda = gobjects(length(Model.dxLambda),1); % Initialize handles
    for idx = 1:length(Model.dxLambda)
        filename = [outputdir '/lambda_' num2str(Model.dxLambda(idx)) 'm.jpeg'];
        hfig_lambda(idx) = showlambda(lambdaResults{idx},timescale1yr,...
            manualcounts,Model,filename);
    end
end

%% Calculate new layer templates:
% Based on data in complete interval.

% In case of floating preprocessing steps, this must be done first, and is
% then based on an average value of lambda over interval.
% Preprocessing specs:
meanLambda = mean(diff(Layerpos.final));
[preprocsteps,preprocdist] = ...
    setpreprocdist(Model.preprocsteps(:,2),meanLambda);
% Preprocess batch data:
data_out = makedatafile(Data.data,Data.depth,preprocsteps,Model.derivatives);

% Derived layer templates:
[Template_New,TemplateInfo_New] = ...
    layerstructure(data_out,Data.depth,Layerpos.final,[],Model,Runtype);

% Plot new layer templates (and compare to original ones):
if Runtype.plotlevel > 0
    color = [1 0 1; 1 1 0; 1 0.5 0.5; 0 0 1]*0.5;
    filename = [outputdir '/layertemplates_new'];
    plotlayertemplates(Template_new,TemplateInfo_new,Model,...
        hfig_template,color,filename);
end

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

%% Calculate similarity to manual layer counts:
% similarity = similarityindex(timescale1yr(:,1),manualcounts(:,1),[],Model.dx,[],Runtype.plotlevel)

%% Show results in matchmaker:
if Runtype.plotlevel>0
    checkinmatchmaker(outputdir,Model,[],Runtype);
end