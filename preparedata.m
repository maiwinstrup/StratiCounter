%% preparedata
% 1) Make a file of manual counts in the correct format
% 2) Prepare data files
% Mai Winstrup, 2014
clear; clc; close all

sett_NEEM_maincore

%% Manual layer counts:
counts = importdata('./Manualcounts/NEEM_main_ALCms.txt');
manualcounts(:,1) = counts.data(:,2);
manualcounts(:,2) = counts.data(:,1); % AD
manualcounts(:,3) = 0; % No uncertainty
manualcounts(:,4) = 0;

%% Save file:
filename = Model.pathManualCounts;
fid = fopen(filename, 'w');

% Add header:
header = ['%% Manual layer counts (' Model.manCountsName ')\r\n\r\n%% Depth '...
    '\t Age \t Uncertainty per layer \t Accumulated uncertainty \r\n'];
fprintf(fid, header);
% Add data:
nDigits = 3;
timescaleformat = ['%.' num2str(nDigits) 'f  \t %.1f \t %.1f \t %.1f']; % First: space, second: #digits

fprintf(fid, [timescaleformat '\r\n'], manualcounts'); 

% Close file:
fclose(fid); 

%% Data files:
data.name = {'NH4','HNO3','BCconc','BCgeom','Na','Mg','nssS', 'ssSNa','Cl',...
    'nssCa','Mn','Sr'};
% mgl: i, Lree

rawdata = xlsread('../Data/NEEM/NEEMSC_database_100314_subset_120514.xlsx');

data.depth = rawdata(:,1);
data.data(:,1) = rawdata(:,5);
data.data(:,2) = rawdata(:,6);
data.data(:,3) = rawdata(:,10);
data.data(:,4) = rawdata(:,11);
data.data(:,5) = rawdata(:,13);
data.data(:,6) = rawdata(:,14);
data.data(:,7) = rawdata(:,17);
data.data(:,8) = rawdata(:,18);
data.data(:,9) = rawdata(:,19);
data.data(:,10) = rawdata(:,21);
data.data(:,11) = rawdata(:,24);
data.data(:,12) = rawdata(:,29);

data.data(data.data==-0.999)=nan;

save ./Data/NEEMdata data
