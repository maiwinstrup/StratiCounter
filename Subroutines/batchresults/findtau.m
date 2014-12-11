function [tau, postau] = findtau(FBprob,meanLambda,d,batchLength,Model,plotlevel)

%% [tau,postau] = findtau(FBprob,meanLambda,d,batchLength,Model,plotlevel)
% Selecting a good pixel for ending the current batch. Successive batches 
% are slightly overlapping with a maximum overlap section depending on the 
% value of Model.batchOverlap (measured in mean layer thicknesses). In this 
% way, any fuzzy annual layering in the last part of the batch is 
% determined also based on subsequent data. 
% tau is selected as the first pixel after a very likely layer boundary, 
% for which we are sure to be out of the previous layer. 
% Furthermore is calculated the probability of where ends the last layer 
% before the one in pixel tau (postau).

% Mai Winstrup
% 2014-10-13 20:03: Cleaned up as separate script
% 2014-10-14 14:56: Plot added, using tnext (previously: tmax+dmax)
% 2015-10-15 16:24: Added calculation of postau
% 2014-10-20 11:47: New version of lastlayerpos

%% Sometimes we wish to use the entire data section: 
% If we have tiepoints, or if batchOverlap is equal to zero. 
if Model.batchOverlap == 0 || ~isempty(Model.tiepoints)
    tau = batchLength; % Using the entire data section
    postau = [];
    return
end

%% Selecting a pixel with no/little probability of being a layer boundary, 
% located right after a very accurately positioned layer boundary.
% Mean, min, and maximum layer thickness in pixels: 
dmean = meanLambda/Model.dx;
dmin = min(d);
dmax = max(d);

% Pixels with a high value of eta_all is indicative of a very accurately 
% positioned layer boundary:
eta_all = sum(FBprob.eta_bar,2);

% Find the pixel with maximum layer ending probability within the selected
% last part of the batch:
tstart = max(1,batchLength-round(Model.batchOverlap*dmean)-dmax);
tend = batchLength-dmax; 
% Find peak locations:
[peakvalue, ipeak] = findpeaks(eta_all(tstart:tend),'minpeakdistance',dmin);
% The maximum peak:
[~, ipeakmax] = max(peakvalue);
tmax = ipeak(ipeakmax)+tstart-1; % [pixel]
% And the subsequent peak: 
if ipeakmax < length(peakvalue)
    tnext = ipeak(ipeakmax+1)+tstart-1; % [pixel]
else
    [~, ipeak] = findpeaks(eta_all(tend+1:end),'minpeakdistance',dmin);
    tnext = ipeak(1)+tend;
end

% We want to pick the first pixel between these two peaks, having (almost) 
% no probability of being a layer boundary:
zerolimit = 10^-3;
tau = find(eta_all(tmax:tnext)<zerolimit,1,'first')+tmax-1;

% Or, if that is not possible, use the best possible location:
if isempty(tau)
    [~, itau] = min(eta_all(tmax:tnext));
    tau = itau+tmax-1; % [pixel]
end

%% Probabilities corresponding to ending location of last layer:
% "Last layer" is the layer previous to the one in pixel tau.
zerolimit = 10^-3;
postau = lastlayerpos(FBprob,tau,d,zerolimit,plotlevel); % 2014-10-20 11:45

%% Plot selected value of tau:
if plotlevel>1
    figure;
    hline(1)=plot(eta_all,'-b');
    hold on
    plot(tmax, eta_all(tmax),'ob')
    plot(tnext, eta_all(tnext),'ob')
    hline(2)=plot(tau,eta_all(tau),'.r','markersize',15);
    hline(3)=plot((tau-length(postau)-1):tau,[0;postau;0],'-r');
    xlim([tstart batchLength])
    xlabel('Pixel')
    ylabel('Probability')
    title('Probability of ending any layer j in given pixel')
    legend(hline,{'sum_j(\eta_{j})','\tau','pos(\tau)'})
end