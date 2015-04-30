function t0 = initialage(manualcounts,layerprobdist,Model)

%% t0 = initialage(manualcounts,layerprobdist,Model)
% Estimate/find age of the first autocounted annual layer boundary, 
% according to the manual counts. While t0 can often be inferred directly, 
% only an estimate can be provided if manual counts do not cover the first 
% layer boundary. The two are then tied at the first manual layer count
% within interval of automated layer counts (if such exists). 
% Value of t0 is in ageUnitOut (as is manualcounts).

% Copyright (C) 2015  Mai Winstrup

%% Age for initial layer boundary:
%% Situation A: Tiepoints are provided.
if ~isempty(Model.tiepoints)
    % Initial age is found from age of first tiepoint:
    t0 = Model.tiepoints(1,2);
    return
end

%% Situation B: Tiepoints are not provided. Initial age is according to 
% the manual counts:
if isempty(manualcounts)
    % If no layer counts: Start with layer 1
    t0=1;
    
elseif manualcounts(1,1)>layerprobdist(end,1)
    % We have no point at which the two timescales overlap, and hence can 
    % be tied.
    t0=1;
    
else
    % The depth in which the two timescales are to be tied:
    dtie = max(manualcounts(1,1)+Model.dx/2,layerprobdist(1,1));

    % Manual age at tiepoint:
    tmanual = floor(interp1(manualcounts(:,1),manualcounts(:,2),dtie));
    
    % Automated layer count number at dtie:
    tauto = interp1(layerprobdist(:,1),layerprobdist(:,2),dtie);
    % No need for rounding after, since all pixels are given a value.

    % Age at start of data series, according to manual counts:
    switch Model.ageUnitOut
        case 'AD'
            t0 = tmanual+tauto;
        otherwise
            t0 = tmanual-tauto;
    end
end