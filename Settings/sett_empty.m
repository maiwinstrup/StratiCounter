%% Example of settings for an ice core record
% Mai Winstrup, 2015

%% Data files:
Model.icecore = 'icecore';
% Data existing for core:
Model.species = {'Isotopes'};
Model.nSpecies = length(Model.species);

% Weighting of various species:
Model.wSpecies = ones(Model.nSpecies,1);

% Path to data file:
Model.path2data = './Data/data_empty.mat';

%% Depth interval [m]:
Model.dstart = 0;
Model.dend = 100; 

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
Model.dx_offset = 0;
% If using e.g. midpoints of dx intervals, the value of dx_offset should
% be set as 0.5 (to avoid unnecessary interpolation). 

% Preprocessing of each data series:
for j = 1:Model.nSpecies
    Model.preprocsteps{j,1} = {'zscore',0.5};
    Model.preprocsteps{j,2} = [];
end
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

%% Length of each data batch (in approximate number of layers):
Model.nLayerBatch = 50; 
% If tiepoints are given, the length of each data batch corresponds to the
% interval between these.

%% Provide path to manual layer counts to be used for initialization:
Model.nameManualCounts = 'counts_empty.txt';
Model.ageUnitManual = 'layers';
% Format of file with manual layer counts:
% counts(:,1): Depth
% counts(:,2): Age
% counts(:,3): Uncertainty of layer (0: certain, 1: uncertain)
% counts(:,4): Accumulated uncertainty from start of interval

%% The annual layer model: 
% Depth interval used for determining layer shapes based on preliminary 
% manual layer counts: 
Model.manualtemplates = [Model.dstart Model.dend];

%% Initial model parameters and variation allowed:
% The initial set of layer parameters will be based on manual layer counts.
% Depth interval used to estimate these:
Model.initialpar = [Model.dstart, Model.dstart+0.5*(Model.dend-Model.dstart)];
% Using the first half of data.

%% Iterations and convergence:
% Number of iterations per batch:
Model.nIter = 4;

% Parameters allowed to be updated at each iteration.
% The ordering is: 
% 1: my, 2: sigma (mode and variace of layer thickness distribution)
% 3: par, 4: cov, 5: nvar (layer shape mean parameters, inter-annual 
% covariance and white noise component)
Model.update = {'ML', 'ML', 'ML', 'ML', 'ML'}; 
% Options:
% 'none': No updates (i.e. maintained as layerpar0)
% 'ML': Maximum-Likelihood updates

%% Output of algorithm:
% Length of interval(s) for determining average layer thicknesses:
Model.dxLambda = [1 5]; % [m]
% If empty, lambda values are not determined.

% Depth sections for calculation of confidence intervals of the number of 
% layers between these marker horizons:
Model.dMarker = [];
% Multiple sets of sections can be included as follows:
% Model.dMarker{1} = [101, 152.5, 204];
% Model.dMarker{2} = [121, 142, 201];
% Distributions are also calculated for the section from the beginning of 
% the data series to the first marker horizon, and from the last
% marker horizon within data series to end of the data series. 

% Which timescale terminology to be used for output? 
Model.ageUnitOut = 'layers';
% Options: AD, BP, b2k, layers