function result = yesnoinputwdefault(promt,default)

coder.extrinsic('input');

%% result = yesnoinputwdefault(promt,default)
% If the input is not n/no/n* or y/yes/y*, it is given the default option (which 
% may either be yes/y or no/n)
% Copyright (C) 2015  Mai Winstrup
% 2014-06-14 13:33

%% Input:
result = input([promt ' (y/n) '],'s');

%% If empty: Use default
if isempty(result); result = default; disp(result); return; end

%% Convert to full/short names for yes/no:
if strcmp(default,'yes') || strcmp(default,'no')
    if strcmp(result(1),'n'); result = 'no'; end
    if strcmp(result(1),'y'); result = 'yes'; end
elseif ismember(default,{'y','n'})
    if strcmp(result(1),'n'); result = 'n'; end
    if strcmp(result(1),'y'); result = 'y'; end
else
    disp('Default value is not recognized')
end

%% For all other input than y/yes/n/no, the default value is used:
switch default
    case 'y'
        if ~strcmp(result,'n'); result = default; end
    case 'yes'
        if ~strcmp(result,'no'); result = default; end
    case 'n'
        if ~strcmp(result,'y'); result = default; end
    case 'no'
        if ~strcmp(result,'yes'); result = default; end
end