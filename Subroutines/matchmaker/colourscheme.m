function [groupspecies, groupcolour] = colourscheme()

%% Define pre-set colour scheme to use in matchmaker datafiles.
% Copyright (C) 2016  Mai Winstrup

%% Common endings and starts that should not affect colour choice:
commonsuffix = {'_DRI','_ic','_deconv','_redo','_1','_2','_3','_repeat','_a','_b','_c'};
commonprefix = {'nss','ss'};

%% Water stable isotopes, d18O: Dark grey
groupspecies{1} = {'d18O','d18O_d','isotopes','Isotopes','d18Odetails','d18O_deconv2'};
groupcolour{1} = [0.1, 0.1, 0.1]; 

% Water stable isotopes, dD: Grey
groupspecies{2} = {'dD','Dexcess'};
groupcolour{2} = [0.3, 0.3, 0.3]; 

%% ECM: Light blue
groupspecies{3} = {'ECM','AC-ECM','DC-ECM'};
groupcolour{3} = [0.5, 0.8, 1];

% DEP and conductivity: Dark blue
groupspecies{4} = {'DEP','DEP adj','cond','Cond'};
groupcolour{4} = [0.3, 0.5, 1];

%% Dust: Light grey
groupspecies{5} = {'Partcnt_m','Partcnt_1','partcnt_m','partcnt','Dust','DustN'};
groupcolour{5} = [0.6, 0.6, 0.6];

% Calcium: Dark green
groupspecies{6} = {'Ca'};
groupcolour{6} = [0.2, 0.5, 0.2];

% Magnesium: Yellowish green
groupspecies{7} = {'Mg'};
groupcolour{7} = [0.7, 0.85, 0.06];

% Rare earth elements: Green
groupspecies{8} = {'Li','Mn','Ce','Pb','LightREE','Sr','LREE','exPb'};
groupcolour{8} = [0.2, 0.7, 0.5];

%% Sodium: Purple
groupspecies{9} = {'Na','sodium'};
groupcolour{9} = [0.86, 0, 1];

% Cloride: Violet
groupspecies{10} = {'Cl'};
groupcolour{10} = [0.6, 0.5, 1];

% Flouride: Bright blue/turquise
groupspecies{11} = {'F'};
groupcolour{11} = [0, 0.86, 0.88];

% Kalium: Blue/purple
groupspecies{12} = {'K'};
groupcolour{12} = [0.3, 0.3, 1];

%% Sulfate-derived species: Orange
groupspecies{13} = {'S','SO4','sulfate','Sulfate'};
groupcolour{13} = [1, 0.5, 0];

% MSA: Light orange
groupspecies{14} = {'MSA'};
groupcolour{14} = [0.9, 0.7, 0.25];

%% Ammonium: Green
groupspecies{15} = {'NH4'};
groupcolour{15} = [0, 1, 0];

%% Nitrate: Red
groupspecies{16} = {'NO3','HNO3'};
groupcolour{16} = [1, 0, 0];

%% Black carbon: Black
groupspecies{17} = {'BC','BCgeom','BCconc','BCgeom3','BCconc3'};
groupcolour{17} = [0, 0, 0];

%% 10Be:
groupspecies{18} = {'Be10','10Be','10Be_average'};
groupcolour{18} = [0.5, 0.2, 0.1];

%% Add suffix and prefix:
for i = 1:length(groupspecies)
    % Add prefix:
    K = length(groupspecies{i});
    for j = 1:length(commonprefix)
        % New names:
        clear newnames
        for k = 1:K
            newnames{k} = [commonprefix{j} groupspecies{i}{k}];
        end
        groupspecies{i} = [groupspecies{i}, newnames];
    end
    % Add suffix:
    K = length(groupspecies{i});
    for j = 1:length(commonsuffix)
        % New names:
        clear newnames
        for k = 1:K
            newnames{k} = [groupspecies{i}{k} commonsuffix{j}];
        end
        groupspecies{i} = [groupspecies{i}, newnames];
    end
end            