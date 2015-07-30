function [probdistSections,layer0new_no,dSectionBounds] = ...
    sectionlayerprobs(depth,FBprob,dSectionBounds,layer0_no,tau,ntau,...
    d,pd,logb,plotlevel)

%% [probdistSections,layer0new_no,dSectionBounds] = sectionlayerprobs(depth,...
%    FBprob,dSectionBounds,layer0_no,tau,ntau,d,pd,logb,plotlevel)
% Calculate mumber of layers and associated uncertainties for sections of
% data. Depth boundaries of the individual sections are given in the 
% variable "dSectionBounds", and sections may both be regular (as used for 
% calculating annual layer thicknesses over equidistant intervals) and 
% irregular (providing layer number distributions for counting between 
% marker horizons).

% Output:
% probdistSections(:,1): Layer numbers
% probdistSections(:,2:M+1): Corresponding probabilities for the M sections
% dSectionBounds: End depth of section bounds ending within current batch (in meters)
% layer0new_no: Layer number probability distribution by end of batch (tau)
% (used for calculating results for next batch)
% Last layer position probabilities used by start of batch are similar to
% those corresponding to full batch. 

% Copyright (C) 2015  Mai Winstrup

%% Location of section boundaries (in pixel) within current batch, up to 
% and including pixel tau:
pxSectionBounds = interp1(depth(1:tau),1:tau,dSectionBounds,'nearest',nan);
% Not using boundaries located in pixel 1, which are taken to be part of 
% the previous batch:
mask = pxSectionBounds>1 & pxSectionBounds<=tau;
pxSectionBounds = pxSectionBounds(mask);

% Corresponding depths in meters:
dSectionBounds = dSectionBounds(mask);

% Number of section boundaries within batch:
nSections = length(pxSectionBounds);

%% In case of no section boundaries within batch:
% The layer number distribution in pixel tau is convolved with the layer 
% number distribution from 'Layer0': 
if isempty(pxSectionBounds)
    % Calculate resulting probabilities for pixel tau, and save these for 
    % initialization of the next batch:
    zerolimit = 10^-4;
    layer0new_no = convolveprobdist(layer0_no,...
        [(0:ntau(end,1))',FBprob.gamma(tau,1:ntau(end,1)+1)'],zerolimit);
    
    % No results for the non-existing section boundaries within interval:
    probdistSections = []; 
    return
end

%% Initializing the array of probability distributions for all sections
% within batch:
% Minimum and maximum number of layers within all sections in batch:
% Individual sections may be smaller/larger than others, and they may start 
% within the batch. Hence, we do not know the lower limit of layers. 
yrmin = 0;
yrmax = layer0_no(end,1)+ntau(end,1); 

% Layer number distributions for all sections in batch:
probdistSections(:,1)=yrmin:yrmax; 
probdistSections(:,2:nSections+1)=0;

%% First section in batch is special: 
% Simply convolve the calculated probability distribution at first boundary 
% with the probabilities corresponding to pixel 1.

% Probabilities for batch at first boundary (pxSectionBounds(1)): 
probdist1 = [(0:ntau(end,1))',FBprob.gamma(pxSectionBounds(1),1:ntau(end,1)+1)'];
% Combine the layer distributions:
zerolimit = 10^-4; 
% Less accuracy is needed than for timescale results: We may use a higher 
% cut-off value.
probdistSection1 = convolveprobdist(layer0_no,probdist1,zerolimit); 

% Inset in layer probability array for all intervals:
probdistSections(probdistSection1(:,1)+1,2) = probdistSection1(:,2);

% Probability corresponding to the location of last layer boundary (i.e. 
% the layer previous to the one in pixel pxSectionBounds(1)):
zerolimit = 10^-3;
layerX_pos = lastlayerpos(FBprob,pxSectionBounds(1),d,zerolimit,[],plotlevel); 

%% Remaining sections in batch:
% Layer thickness parameters:
dmax = max(d);
D = length(d);
% Total length of batch:
batchLength = size(FBprob.gamma,1);

for iSection = 2:nSections+1 % Also for the last section, that runs until end of batch
    %% Start pixel of interval:
    t0 = pxSectionBounds(iSection-1);
 
    %% Form the initial condition vector, layer0section_pos:
    % Ending location of previous layer relative to t0 is:
    layer0section_pos = zeros(dmax,1);
    layer0section_pos(end-length(layerX_pos)+1:end)=layerX_pos;
    
    %% If 2 options: selecting the most likely location of the layer
    % boundary?
    
    %% Obtain layering (and associated uncertainties) within section:
    % Length of section:
    % Using all rest of batch, so that all subsequent data is used
    sectionLength = batchLength-t0+1; 
    % Alternatively, could also use: 
    % sectionLength = t1-t0+1; 
    
    % Maximum number of layers in section:
    nMinStart = find(FBprob.gamma(t0,:)>0,1,'first');
    nMaxEnd = find(FBprob.gamma(end,:)>0,1,'last');
    % If using other definition of sectionLength, then here use: 
    % nMaxEnd = find(FBprob.gamma(t1,:)>0,1,'last');
    nLayerMaxSection = nMaxEnd-nMinStart+1;
    
    % Run Forward-Backward with the above as input:
    [~,~,~,~,gammasection,~,~,~]=...
        forwardbackward(logb(t0:end,:),log(pd),layer0section_pos,...
        sectionLength,d,dmax,D,nLayerMaxSection,[],plotlevel);    
 
    %% Results for section: 
    if iSection<=nSections
        % End pixel of interval: Both endpoints are included in interval
        t1 = pxSectionBounds(iSection);
        
        % Use the calculated gamma distribution at t1: 
        clear layerprob
        layerprob(:,1) = 0:nLayerMaxSection-1; 
        layerprob(:,2) = gammasection(t1-t0+1,:);
        % Compact distribution: 
        zerolimit = 10^-3;
        layerprob = compactprobdist(layerprob,zerolimit);
        
        % Inset in layer probability array for all intervals:
        probdistSections(layerprob(:,1)+1,1+iSection) = layerprob(:,2);

        if t1 == tau
            % This is the last section in batch; next section starts in 
            % pixel tau, which is the same as pixel 1 in next batch. 
            % Next section in next batch restarts with layer 0:
            layer0new_no(1,:) = [0, 1];
            break
        else
            % New layerX_pos to be used as input for subsequent section:
            zerolimit = 10^-3;
            layerX_pos = lastlayerpos(FBprob,t1,d,zerolimit,[],plotlevel); 
        end
        
    elseif iSection==nSections+1
        % We are in the very last section, which is only partly covered by
        % the current batch. 
        
        % Resulting probability distribution at pixel tau is used as
        % initial condition for next batch:
        % Calculate the gamma distribution at tau, to be used as input 
        % for layer0 for next batch: 
        ntausection(:,1) = 0:nLayerMaxSection-1; 
        ntausection(:,2) = gammasection(tau-t0+1,:);
        
        % Remove tails of distribution etc:
        zerolimit = 10^-3;
        ntausection = compactprobdist(ntausection,zerolimit);
       
        % Inset as initial condition for next batch:
        layer0new_no = ntausection; 
    end
end
    
%% Remove zero entries:
index1 = find(sum(probdistSections(:,2:end),2)>0,1,'first');
index2 = find(sum(probdistSections(:,2:end),2)>0,1,'last');
probdistSections = probdistSections(index1:index2,:);