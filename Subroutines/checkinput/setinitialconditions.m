function [batchStart,Layer0,Template,Prior,Layerpar,dDxLambda,...
    logPobs,logPobsNorm,relweight,Result] = setinitialconditions(Data,...
    Model,counts0,Template0,Layerpar0,nBatch)

%% [batchStart,Layer0,Template,Prior,Layerpar,dDxLambda,logPobs,...
%    logPobsNorm,relweight,Result] = setinitialconditions(Data,Model,...
%    counts0,Template0,Layerpar0,nBatch)
% Initialize matrices and set values corresponding to initial conditions.
% Matrices for results are also initialized. 
% Copyright (C) 2015  Mai Winstrup

%% Initialize empty matrices:
[batchStart,Layer0,Template,Prior,Layerpar,logPobs,logPobsNorm,...
    relweight,Result] = initializematrices(nBatch,Model);

%% Start pixel for first batch: 
if ~isempty(Model.tiepoints)
    % Start in pixel corresponding to the uppermost tiepoint:
    batchStart(1) = interp1(Data.depth,1:length(Data.depth),Model.tiepoints(1,1),'nearest'); 

elseif isempty(counts0)||counts0(1,1)>Model.dstart+1;
    % If manual counts are not provided (for the entire interval, or for 
    % the first 1 meter of the data).
    batchStart(1) = 1;
else
    % Depth of the first (certain) layer in the manually counted chronology:
    firstcertainlayer = find(counts0(:,3)==0,1,'first');
    dfirstcertainlayer = counts0(firstcertainlayer,1); % [m]
    % Pixel corresponding to the start of a new layer:
    batchStart(1) = interp1(Data.depth,1:length(Data.depth),dfirstcertainlayer+Model.dx/2,'nearest'); 
end

%% Probabilities for location and layer number of first layer boundary:
% Location of the first layer boundary is taken to have negligible
% uncertainty, i.e. its associated probability distribution is given by:
Layer0(1).pos = 1; 

% Probability of layer number for the previous layer (that just ended):
Layer0(1).no(:,1)=0; % Layer number
Layer0(1).no(:,2)=1; % Probability

%% Intervals for calculation of average layer thicknesses:
% Number of different partitions of the data series to be used for mean 
% layer thickness calculations:
nDxLambda = length(Model.Out.dxLambda);

% Boundaries associated with the depth sections:
dDxLambda = cell(1,nDxLambda);
for i = 1:nDxLambda
    % Interval boundaries:
    d0 = floor(Data.depth(batchStart(1))/Model.Out.dxLambda(i))*Model.Out.dxLambda(i); 
    dDxLambda{i} = d0:Model.Out.dxLambda(i):Data.depth(end);   
    % First value may be outside interval. If so, it is removed:
    dDxLambda{i} = dDxLambda{i}(dDxLambda{i}>=Data.depth(batchStart(1)));
end

% Same initialization as for first layer:
for i = 1:nDxLambda
    Layer0(1).noDx{i} = Layer0(1).no;
end

%% Marker horizons:
% Same initialization as for first layer:
for i = 1:length(Model.Out.dMarker)
    Layer0(1).noMarker{i} = Layer0(1).no;
end

%% Initialize layer templates and layer parameters:
% Initial templates for batch 1:
Template(:,1,1) = Template0;

%% Prior on parameter estimates:
% The values in Layerpar0 are used as a priori estimates for the first
% batch of data, and are hence our initial guess of the layer parameters. 

% Layer thicknesses: 
% Mode of distribution:
Prior(1,1).m = Layerpar0.my;
% A priori knowledge on mode value (only if using QB-reestimation):
if strcmp(Model.update{1},'QB')
    % Theoretical uncertainty of the mean, as if based on a single previous
    % batch:
    % theoretical_var_of_mean = Layerpar0.sigma^2/Model.nLayerBatch;
    Prior(1,1).v = 0.5; 
    disp(['Using initial prior: Prior.v = ' num2str(Prior(1,1).v)])
end
% Interannual variation in layer thicknesses:
Prior(1,1).sigma = Layerpar0.sigma;

% Layer parameters:
Prior(1,1).u = Layerpar0.par;
% A priori knowledge on mode value (only if using QB-reestimation):
if strcmp(Model.update{3},'QB')
    % Theoretical uncertainty of the template parameter, as if based on a
    % single previous batch:
    % theoretical_cov_of_mean = Layerpar0.cov/Model.nLayerBatch;
    Prior(1,1).invU = inv(1);
    disp(['Using Prior.invU = ' num2str(Prior(1,1).invU)])
end
% Interannual variation in layer parameters:
Prior(1,1).cov = Layerpar0.cov;
% White noise component:
Prior(1,1).nvar = Layerpar0.nvar;

% For gibbs sampling:
if sum(strcmp(Model.update,'gibbs')>0)
    % Prior for my: variance
    Prior(1,1).kappa = 1/0.1^2; 

    % Prior for sigma: Gamma distribution
    Prior(1,1).alpha = 10;
    Prior(1,1).beta = 0.4/1.5;
    
    % Precision matrix for parameter vector:    
    for j = 1:Model.nSpecies
        for i = 1:Model.order
            Prior(1,1).Lambda(i,i,j) = abs(((Prior(1,1).u(i,j))/2)^-1); % very loose prior
        end
    end
end

%% Layer parameters: Initial values for batch 1:
Layerpar(1,1,1) = Layerpar0;

%% Initialize matrices for performance evaluation: 
% Evolution of log(Pobs) with iterations:
logPobs(1,1,1) = -inf;