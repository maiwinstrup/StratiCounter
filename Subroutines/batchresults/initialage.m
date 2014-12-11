function t0 = initialage(manualcounts,layerdist,Model)

%% t0 = initialage(manualcounts,layerdist,Model)
% Estimate/find age of the first autocounted annual layer boundary, 
% according to the manual counts. While t0 can often be inferred directly, 
% only an estimate can be provided if manual counts do not cover the first 
% layer boundary. The two are then tied at the first manual layer count
% within interval of automated layer counts (if such exists). 
% Value of t0 is in ageUnitOut (as is manualcounts).

% Mai Winstrup
% 2014-10-21 15:49: First version

%% Age for initial layer boundary, according to manual counts:
if isempty(manualcounts)
    % If no layer counts: Start with layer 1
    t0=1;
    
elseif manualcounts(1,1)>layerdist(end,1)
    % We have no point at which the two timescales overlap, and hence can 
    % be tied.
    t0=1;
    
else
    % The depth in which the two timescales are to be tied:
    dtie = max(manualcounts(1,1)+Model.dx/2,layerdist(1,1));

    % Manual age at tiepoint:
    tmanual = floor(interp1(manualcounts(:,1),manualcounts(:,2),dtie));
    
    % Automated layer count number at dtie:
    tauto = interp1(layerdist(:,1),layerdist(:,2),dtie);
    % No need for rounding after, since all pixels are given a value.

    % Age at start of data series, according to manual counts:
    switch Model.ageUnitOut
        case 'AD'
            t0 = tmanual+tauto;
        otherwise
            t0 = tmanual-tauto;
    end
end        