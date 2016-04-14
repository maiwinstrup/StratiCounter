%% prepareinputfiles
% This script will help you to prepare input files for StratiCounter. 
% The script will have to be modified according to the format of the input 
% data.
% Mai Winstrup, 2016
clear; clc; close all

%% Construct the data file structure array:
% If importing data from excel sheet: 
rawdata = xlsread('path2datasheet.xls'); % replace with proper path to file

% Or: Import data from txt-file:
rawdata = importdata('path2datasheet.txt'); % replace with proper path to file

% Now you can build a structure array containing the four fields: name, 
% data, depth, and depth_no, as required by the StratiCounter script. 
% This is done using dot-notation (structureName.fieldname).  

% Assume the first column of data file is depth. 
% In this case, the first column is assigned to the field "depth". 
data.depth{1} = rawdata(:,1);
% We will in the following assume that all data in the data file 
% corresponds to this depth scale.

% The data records on this depth scale are assinged to the field "data".
% This is done for each of the N data series in file.
for n = 1:N % N records  
    data.data{n} = rawdata(:,n+1); 
    % We here assume that the data starts in column 2, as the first column 
    % contained the depth scale. 
end

% Provide the name of the data series:
data.name{1} = 'name of 1st record';
data.name{2} = 'name of 2nd record';
...
data.name{N} = 'name of Nth record';

% In this simple example, all data records share a common depth scale. This
% depth scale is in data.depth{1}, i.e. the first cell of the depth field. 
% We indicate this relationship by setting all N values in the field 
% "depth_no" equal to 1:
data.depth_no(1:N) = 1;

% You have now produced the data file! 
% Save your datafile in the "Data" folder:
save ./Data/yourdata data 

% In the more complex case where data records are on multiple depth
% scales, these depth scales should be added in additional cells under the 
% field "depth", i.e. data.depth{2} = ..., data.depth{3}=....
% The numbers in the field "depth_no" should then be changed to reflect
% which depth scale corresponding to the N data records. 
% In this case, data will likely also be combined from various data files, 
% each being imported separately. 

%% Produce a text file containing the manual layer picks:
% Import the depths of manual layer picks from either xls-file or txt-file: 
counts = xlsread('path2manualcounts.xls'); % replace with proper path to file
counts = importdata('path2manualcounts.txt'); 

% Include information on the manual layer counts in the following order: 
% First column contains the depth of the manual layer counts:
manualcounts(:,1) = counts(:,1);
% Second column contains the corresponding (estimated) age for the manual
% layer counts. Age may also be given as "layer number". 
manualcounts(:,2) = counts(:,4);
% In the above, ages are assumed to be given in column 4 of the original 
% file with layer counts. 

% The 3rd and 4th columns concerns the estimated uncertainty of the manual
% layer counts, if such uncertainty estimates exist. 
% In the 3rd column, the uncertainty of the individual layer counts is 
% provided. If a layer is deemed uncertain, the corresponding "uncertainty" 
% value should be equal to 1. If not, it is given a value of 0. 
% The 4th column contains the accumulated uncertainty estimates, also where 
% the individual layer uncertainty estimates do not exist. 

% If uncertainty estimates are non-existing, the uncertainty values are set 
% equal to 0:
manualcounts(:,3) = 0;
manualcounts(:,4) = 0;

% You have now produced the file with a set of manual layer counts! 

% This file must now be saved in the folder "Manualcounts":
% Below is shown an example of how it can be saved as textfile with header.
filename = './Manualcounts/yourcorename.txt';
fid = fopen(filename, 'w');

% Add header:
header = ['%% Manual layer counts (' manCountsName ')\r\n\r\n%% Depth '...
    '\t Age \t Uncertainty per layer \t Accumulated uncertainty \r\n'];
fprintf(fid, header);
% Add data:
nDigits = 3;
timescaleformat = ['%.' num2str(nDigits) 'f  \t %.1f \t %.1f \t %.1f']; % First: space, second: #digits
% Save data:
fprintf(fid, [timescaleformat '\r\n'], manualcounts'); 
% Close file:
fclose(fid); 