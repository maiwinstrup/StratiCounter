function Model = defaultsettings()

%% Model = defaultsettings()
% Provide default settings for StratiCounter.
% Copyright (C) 2015  Mai Winstrup

%% Data files:
Model.icecore = '';
% If using synthetic data, the name should be set as 'SyntheticData'.
% Data existing for core:
Model.species = '';
Model.nSpecies = 0; 

% Weighting of various species:
Model.wSpecies = ones(Model.nSpecies,1);

% Path to data file:
Model.path2data = '';

%% Depth interval [m]:
Model.dstart = [];
Model.dend = []; 

%% Marker horizons to be used as tiepoints?
% No tiepoints:
Model.tiepoints = [];
% With tiepoints:
% Model.tiepoints(:,1) = []; % Depth [m]
% Model.tiepoints(:,2) = []; % Corresponding age
% Model.ageUnitTiepoints = ''; % Age unit of tiepoints
% Options: AD, BP, b2k, layers

%% Data treatment:
% Resolution of data series to be used:
Model.dx = 10^-3; % [m/px]
% If Model.dx is empty: No interpolation to equidistant depthscale (this 
% option is not fully implemented). 
% The resolution is allowed to (slowly) change with depth. 
Model.dx_offset = 0;
% If using e.g. midpoints of dx intervals, the value of dx_offset should
% be set as 0.5 (to avoid unnecessary interpolation). 

% Preprocessing of each data series:
Model.preprocsteps = cell(Model.nSpecies,2); 
% 1st row: Initial preprocessing
% 2nd row: Batch preprocessing 

% Possible preprocessing steps:
% Non-distance dependent transformations:
% - No processing ([])
% - Interpolation over nans ('interpNaNs')
% - Subtract mean ('minusmean')
% - Log-transformation ('log')
% - Box-Cox transformation ('boxcox',[],[lambda, (alpha)])
% Distance-dependent transformations:
% - Normalization using standard deviation ('zscore',Lwindow)
% - Normalization using min-max ('minmax', Lwindow)
% - Normalization using quantiles ('quantile', Lwindow)
% - Using CDF-transform ('cdftransform',Lwindow)
% - Subtract baseline ('minusbaseline',Lwindow,(quantile))
% - Subtract smooth curve calculated using running average ('minussmooth', Lwindow)
% Window lengths are given in m (1st row) or in terms of lambda (2nd row). 

% Using data series itself and (possibly) its derivatives: 
% Number to derivative data series to include (apart from original data set)
Model.derivatives.nDeriv = 1; % Use data and first derivatives
% Calculation of slope and curvature using Savitsky-Golay:
% First entry is for 1. derivatives, second entry is for 2. derivatives
Model.derivatives.slopeorder = [1 2]; 
Model.derivatives.slopedist = [3 5]; % Window size (in number of observations)

%% Length of each data batch (in approximate number of layers):
Model.nLayerBatch = 50; 
% If tiepoints are given, the length of each data batch corresponds to the
% interval between these. 

%% If using manual counts for initialization:
Model.nameManualCounts = '';
Model.ageUnitManual = '';
% Format of file with manual layer counts:
% counts(:,1): Depth
% counts(:,2): Age
% counts(:,3): Uncertainty of layer (0: certain, 1: uncertain)
% counts(:,4): Accumulated uncertainty from start of interval

%% The annual layer model: 
% Depth interval used for determining the annual layer model based on 
% preliminary manual layer counts:
Model.manualtemplates = [];
% Not implemented: "[]" should mean to not use manual counts for shape 
% derivations, but to use e.g. a sinusoidal shape.

% Calculation of the emission probabilities (b).
Model.type = 'PCA';
% Not implemented: Other options to be made available (e.g. 'FFT')
Model.order = 1;
% Using polynomial approximations to principal components of order:
Model.pcPolOrder = 5;

% Normalize each layer individually before calculating probabilities?
Model.normalizelayer = 'minusmean'; 
% Options: 
% 'none', 'minmax', 'zscore', 'minusmean'

% Method for calculating the probabilities:
Model.bcalc = 'BLR'; % Bayesian Linear Regression
% Not implemented: Possible options to be made available:
% 'BLR', 'BLR_wNaN', 'physical', 'hierachy', 'fourier', etc.

% Timestep used for up/downsampling layers to stack:
Model.dtstack = 1/64; %[year]

% Number of layer shape iterations:
% Per batch:
Model.nTemplateBatch = 1; 
% For complete data set:
Model.nTemplateFull = 1;
% If 1: Layer shapes are based on manual counts
% If > 1: Layer shapes are based on autocounts from previous counting. 

%% Calculation of p(d):
% Type of layer thickness distribution:
Model.durationDist = 'logn';
% Definition and treatment of tails of distribution:
Model.tailType = 'evenly'; 
% Options: 
% 'Evenly' (remaining tails are evenly distributed; probably best option)
% 'Spliced' (remaining tails spliced onto the respective ends)
Model.tailPrc = 0.5; % Percent removed from each tail. 

% Weighting of b relative to p(d):
Model.bweight = 1;
% Number of species and derivatives are accounted for automatically. 
% A value of Model.bweight>1: The dependency on layer shapes are 
% enhanced relative to layer durations. 

%% Initial model parameters and variation allowed:
% The initial set of layer parameters will be based on manual layer counts.
% Depth interval used to estimate these:
Model.initialpar = []; 

% Using the entire parameter covariance matrix?
Model.covariance = 'none'; 
% Options: 
% 'none': All layer Model parameters are considered uncorrelated.
% 'species': Layer parameters for a specific species is considered 
% correlated, but parameters for individual species are uncorrelated.
% 'full': Both inter- and intra-species correlation is taken into account.

% Entries in noise weighting matrix (W), giving the relative white noise 
% level of derivative data series:
Model.derivnoise = 'analytical';
% Options: 
% 'analytical': Analytical values are used. 
% 'manual': Values based on the manual counting are used. 
% The value of derivnoise can also be ascribed specific numbers. 

%% Iterations and convergence:
% Maximum number of iterations per batch:
Model.nIter = 4;

% Limit for convergence of EM-algorithm: 
Model.eps = -1; 
% If negative, the model will run exactly nIter iterations

% Parameters allowed to be updated at each iteration:
% The ordering is: 
% 1: my, 2: sigma (mode and variace of layer thickness distribution)
% 3: par, 4: cov, 5: nvar (layer shape mean parameters, inter-annual 
% covariance and white noise component)
Model.update = {'ML', 'ML', 'ML', 'ML', 'ML'}; 
% Options:
% 'none': No updates (i.e. maintained as initialpar)
% 'ML': Maximum-Likelihood updates
% 'QB': Quasi-Bayes updates

% Forgetting parameter in QB updates:
if sum(ismember(Model.update,'QB')>0)
    Model.rho = 0.8;
end

%% Combining batches: 
% Interval before end of batch in which the last layer boundary is found:
Model.batchOverlap = 5; % Measured in units of mean layer thicknesses

%% Output of algorithm:
% Percentiles for confidence intervals:
Model.confInterval = [50 95]; 

% Interval(s) for determining average layer thicknesses:
% Regular length intervals:
Model.dxLambda = []; % [m]
% If empty, lambda values are not determined.

% Depth sections for calculation of confidence intervals of the number of 
% layers between these marker horizons:
Model.dMarker = [];

% Which timescale terminology to be used for output? 
Model.ageUnitOut = 'layers';
% Options: AD, BP, b2k, layers