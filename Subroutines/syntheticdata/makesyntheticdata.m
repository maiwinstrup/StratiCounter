 function [Data,counts,Model] = makesyntheticdata(Model,Runtype,outputdir)

%% [Data,counts,Model] = makesyntheticdata(Model,Runtype,outputdir)
% Construct synthetic dataseries according to selected model settings.
% Layer counts are located midway between data points belonging to 
% different layers.

% Mai Winstrup
% 2014-06-20 23:38: Major updates
% 2014-08-19 22:05: shape->template, changes to e.g. outputdir, distinction 
%                   between where layers start and layer counts.
% 2014-08-21 12:42: Outputdir moved to other script, plot templates moved to
%                   checksyntheticdata.m
% 2014-08-21 20:09: Added info on "manual" counts

%% Depth scale: 
Data.depth = (Model.dstart:Model.dx:Model.dend)';
Model.dend = Data.depth(end);
% Length of data series [px]:
dataLength = length(Data.depth);

%% Construct array of true layer boundaries ("counts"):
% Estimate maximum number of years in data series: 
meanLambda = exp(Model.SynthData.Modelpar.my+Model.SynthData.Modelpar.sigma^2/2);
Nmax = 2*round((Model.dend-Model.dstart)/meanLambda);
% Duration of individual layers:
dur = exp((Model.SynthData.Modelpar.my+Model.SynthData.Modelpar.sigma*randn(1,Nmax))); 
% Location of layer starts (depth of first data point in layer):
dlayerstart = Model.dstart+[0 cumsum(dur)];
% Remove layer starts after Model.dend:
dlayerstart(dlayerstart>Model.dend)=[];

% Round layer boundaries to nearest data point:
dlayerstart_px = interp1(Data.depth,1:length(Data.depth),dlayerstart,'nearest');

% In proper units, and with the layer boundary located between data points 
% belonging to different layers: 
counts(:,1) = Data.depth(dlayerstart_px)-Model.dx/2;

% Layers are numbered downwards from the top: 
counts(:,2) = 1:length(counts(:,1));

% No associated uncertainties:
counts(:,3) = 0;
counts(:,4) = 0;
        
%% Duration of layers in pixels:
% Number of layers (including the last fraction of a layer):
nLayer = size(counts,1);

% Layer durations in pixels, corresponding to rounded layer positions:
dur_px = diff(dlayerstart_px);
% Duration of last layer (in pixels):
dur_px(nLayer) = round(dur(nLayer)/Model.dx);

%% Generate parameter values from normal distribution with specified mean 
% vector and covariance matrix:
R = chol(Model.SynthData.Modelpar.cov);
par = repmat(Model.SynthData.Modelpar.par(:)',nLayer,1)+...
    randn(nLayer,Model.nSpecies*Model.SynthData.order)*R;
par = reshape(par,[nLayer Model.SynthData.order Model.nSpecies]);

%% Form data series with successive annual layer signals:
target = nan(dataLength,Model.nSpecies);

% Same number of species:
Model.SynthData.nSpecies = Model.nSpecies;
% Only calculate actual signal, not derivatives:
Model.SynthData.deriv = 0;

% Form layers:
for iLayer = 1:nLayer
    % Construct design matrix for trajectories:
    Z = designmatrix(Model.SynthData,Model.SynthData.Template,dur_px(iLayer)); % 2014-08-18 14:50
        
    % Target section of data series: 
    istart = dlayerstart_px(iLayer);
    % Length of data section covered by layer:
    L = min(dur_px(iLayer),dataLength-dlayerstart_px(iLayer)+1);
    iend = istart+L-1;
    
    % Inset layer signals:
    for j = 1:Model.nSpecies
        % Calculate mean signal:
        x = 1/(2*dur_px(iLayer)):1/dur_px(iLayer):1;
        meansignal = polyval(Model.SynthData.Template(j).mean,x);
        % Target signal is the mean signal plus variation of trajectories:
        target(istart:iend,j) = meansignal(1:L)' + Z(1:L,:,j)*par(iLayer,:,j)';
    end
end

%% Add gaussian white noise:
noise = nan(dataLength,Model.nSpecies);
for j = 1:Model.nSpecies
    noise(:,j) = Model.SynthData.Modelpar.nvar(j)^0.5*randn(1,dataLength);
end

% This noise is added to the target data series, thereby producing a noisy 
% synthetic data sequence:
Data.data(:,1,:) = target + noise;

%% Smooth the resulting data series:
% No smoothing is performed if Model.SynthData.smoothing is empty.
if ~isempty(Model.SynthData.smoothing)
    % Smoothing distance in pixels:
    dx_smoothing = round(Model.SynthData.smoothing/Model.dx);
    for i = 1:Model.nSpecies
        Data.data(:,1,i) = smooth(Data.data(:,1,i),dx_smoothing);
    end
end

%% Plot data and parameter distributions:
if Runtype.plotlevel>0
    Model.SynthData.dx = Model.dx;
    filename = [outputdir '/' 'Data.jpeg'];
    plotsyntheticdata(Data.data,Data.depth,target,par,dur_px,counts(:,1),...
        Model.SynthData,filename,Runtype.plotlevel); % 2014-08-19 20:56
end

%% Calculate derivatives:
for j = 1:Model.nSpecies
    [slope, dslope, wWhiteNoise] = calculateslope(Data.data(:,1,j),Model.slopeorder,Model.slopedist,0); % 2014-07-16 14:06
    % Add to data array:
    Data.data(:,2,j)=slope; 
    Data.data(:,3,j)=dslope;
end

% Add analytical weighting values to model array:
if strcmp(Model.wWhiteNoise,'analytical')
    Model.wWhiteNoise = wWhiteNoise;
end

%% Save data:
save([outputdir '/Data.mat'],'Data','target','par')
save([outputdir '/Layers.txt'],'counts','-ascii')

%% Add info to Model array:
Model.manCountsName = 'Randomly produced synthetic layers';
% Path to the "manual" layercounts.
Model.manCountsPath = [outputdir '/Layers.txt'];
Model.ageUnitManual = 'layers';