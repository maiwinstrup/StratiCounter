function SyntheticModel = checksyntheticmodel(Model,outputdir,plotlevel)

%% SyntheticModel = checksyntheticmodel(Model,outputdir,plotlevel)
% Check for inconsistencies in structure of Model.SynthData, and convert 
% provided mean signal and trajectories to their polynomial approximations.

% Mai Winstrup
% 2014-08-21 14:15: Initial version
% 2014-08-23 12:54: Changes to colors in plot

%% Check for inconsistencies in format of the synthetic "Modelpar" array:
SyntheticModel = Model.SynthData;
if sum(size(SyntheticModel.Modelpar.par)~=[SyntheticModel.order, Model.nSpecies])>0;
    disp('Incorrect format for Model.SynthData.Modelpar.par')
end

if sum(size(SyntheticModel.Modelpar.nvar)~=[SyntheticModel.order, Model.nSpecies])>0;
    disp('Incorrect format for Model.SynthData.Modelpar.nvar')
end

if sum(size(SyntheticModel.Modelpar.cov)~=[SyntheticModel.order*Model.nSpecies,...
        SyntheticModel.order*Model.nSpecies])>0;
    disp('Incorrect format for Model.SynthData.Modelpar.cov')
end

%% Polynomial approximation of mean signal and trajectories:
Template_new(1:Model.nSpecies) = struct('mean',[],'dmean',[],'d2mean',[],...
    'traj',[],'dtraj',[],'d2traj',[]);
for j = 1:Model.nSpecies
    Template_new(j) = polyapprox(SyntheticModel.Template(j).mean,...
        SyntheticModel.Template(j).traj,SyntheticModel); % 2014-08-21 14:08
end
% Replace in SyntheticModel:
SyntheticModel.Template = Template_new;

%% Plot the templates based on which the synthetic data is formed:
if plotlevel>0
    % Add info to SyntheticModel:
    SyntheticModel.icecore = Model.icecore;
    SyntheticModel.species = Model.species;
    SyntheticModel.nSpecies = Model.nSpecies;

    % Plot templates:
    filename = [outputdir '/SyntheticTemplate.jpeg'];
    % Colors for plotting templates for mean (first) and PCs:
    color = [0.5 0.5 0.5; 0 0 1; 0 1 0; 1 0 0];
    plotlayertemplates(SyntheticModel.Template,[],SyntheticModel,[],color,filename); % 2014-08-21 12:56
    
    % Maybe close figure: 
    if plotlevel==1; close(gcf); end
end