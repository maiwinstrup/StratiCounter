function [layerpos,layerpos_issues] = findlayerposfrommode(modevalue)

%% [layerpos, layerpos_issues] = findlayerposfrommode(modevalue)
% Finding the maximum likelihood layer positions based on the
% Forward-Backward algorithm layer-number distributions. 
% Given that the layer probability distributions are based on the MAP 
% criterion, we may have a jump of multiple years occuring in a pixel, and 
% similarly it may happen that we go back to a previous layer (the two 
% events often co-occur). In these cases, we try to resolve the issue as
% "intuitively" as possible; layerpos will contain the best set of 
% boundaries (all corresponding to 1 yr), and layer-positions with issues
% are provided separately thus allowing for manual check-up. 
% Layer positions are given as located midway between datapoints belonging 
% to the two layers. 

% Copyright (C) 2015  Mai Winstrup

%% Find locations with a change in ML layer number:
% Jumps occurring in modevalue:
existingjumps = unique(diff(modevalue(:)))';
existingjumps = existingjumps(existingjumps~=0);
if isempty(existingjumps); 
    layerpos = []; layerpos_issues = [];
    return; 
end

layerpos = [];
jumpvalue = [];
for jump = existingjumps
    pixels = find(diff(modevalue)==jump); % where layers end
    layerpos = [layerpos; pixels(:)];
    jumpvalue = [jumpvalue; jump*ones(length(pixels),1)];
end

% Sort arrays according to depth:
if ~isempty(layerpos)
    [layerpos, index] = sortrows(layerpos,1);
    jumpvalue = jumpvalue(index);
end

% The layer boundary is taken as located midway between data corresponding 
% to the two layers (i.e. is located between two pixels):
layerpos = layerpos+1/2;

%% Do problematic boundaries exist?
% The problematic jumps where the jumpvalue is not equal to +/-1 (value 
% depends on the applied age unit). 
normaljump = mode(jumpvalue);
index = find(jumpvalue'~=normaljump);

% Seemingly "correct" jumps are also problematic if they occur just before 
% or right after a jump of opposite sign to normaljump. We select the one(s)
% closest to any negative jump (if such exists). 
negativejumps = find(sign(jumpvalue)~=normaljump);

if ~isempty(negativejumps)
    % Which jump is closest to the negative jump?
    % Extending layerposition to avoid problems with margins:
    layerpos_ext = [-inf; layerpos; inf];
    % (Adding +1 to all due to the extension above)
    dist_before = layerpos_ext(negativejumps+1)-layerpos_ext(negativejumps);
    dist_after = layerpos_ext(negativejumps+2)-layerpos_ext(negativejumps+1);

    % Selecting the one having the smallest distance:
    for i = 1:length(negativejumps)
        if dist_before(i) >= dist_after(i)
            index = [index, negativejumps(i)+1];
        else
            index = [index, negativejumps(i)-1];   
        end
    end
    
    % Sorting, and ensuring we are within limits of data section:
    index = sort(index);
    index = index(index>=1|index<=length(layerpos));
    % Remove values occuring more than once:
    index = unique(index);
end

layerpos_issues = [layerpos(index), jumpvalue(index)];

%% The following is only performed if problematic boundaries do occur 
% (which is rarely the case).
if isempty(layerpos_issues); return; end

%% Comparing the troublesome boundaries in sections:
% If only one boundary exists:
if length(index) == 1
    % If opposite sign of normaljump, the boundary is removed (this should 
    % never happen, since we would then have added an extra surrounding
    % layer boundary above):
    if sign(jumpvalue(index))~=sign(normaljump)
        layerpos(index) = nan;
    end
    % If of same sign, the boundary is kept. (In this case, the boundary
    % actually ought to be counted as 2 or more layers, but we don't know 
    % where to add the extra layer(s)).

else
    % If more than one boundary exist:
    % How many problematic sections do we have, and where in "index" do the 
    % sections start and end?
    istart = [1 find(diff(index)>1)+1];
    iend = [istart(2:end)-1, length(index)];
    % Converting to index number:
    istart = index(istart);
    iend = index(iend);
    % Number of problematic sections:
    nsections = length(istart);

    for i = 1:nsections
        indexsection = istart(i):iend(i);
    
        % How many layers do the boundaries add up to? 
        jumpsum = sum(jumpvalue(indexsection));
        if jumpsum == 0
            % All layer boundaries are removed:
            layerpos(indexsection) = nan;
    
        elseif sign(jumpsum) ~= sign(normaljump)
            % In total, they add up to a negative number of layers, and 
            % thus, they should all be removed:
            layerpos(indexsection) = nan;
            % This should not happen often, but would in this case cause us 
            % to have added an extra layer boundary somewhere. 

        elseif jumpsum == normaljump
            % Layer boundaries with opposite sign of "normaljump" are removed:
            mask = sign(jumpvalue(indexsection))~=sign(normaljump);
            jumpvalue(indexsection(mask)) = nan;
        
            % The layer boundary with the largest jump value (only those of 
            % the right sign are left) is kept, all others are removed:
            [~,imax] = max(abs(jumpvalue(indexsection))); 
            % If the maximum value occurs more than once, imax only gives 
            % its first entry. Hence we are certain to end up with a single 
            % layer boundary, just as we wish.
            mask=true(length(indexsection),1);
            mask(imax)=false;
            jumpvalue(indexsection(mask)) = nan;
        
        else
            % In this case, either: 
            % A) We have a single layer boundary with a jump value of more 
            % than 1. In this case, we wish to keep it (we will then in 
            % total have missed a layer boundary somewhere), or 
            % B) We have two bad sections are located close to each other, 
            % and is hence regarded as a combined area. The boundaries need 
            % to be combined into a number of layers equal to "jumpsum".
        
            % Situation A: Keep boundary and proceeed to next section. 
            if length(indexsection)==1
                continue
            end

            % Situation B:
            % The layer boundaries with the largest jump values are kept, 
            % all others are removed:
            [~,ilargestjumps] = sort(sign(normaljump)*jumpvalue(indexsection),1,'descend');
            % Keeping only the appropriate number of largest jumps:
            ilargestjumps = ilargestjumps(1:abs(jumpsum));
        
            mask=true(length(indexsection),1);
            mask(ilargestjumps)=false;
            jumpvalue(indexsection(mask)) = nan;
        end
    end
end

%% Remove NaNs in layerpos:
layerpos = layerpos(isfinite(layerpos));