function lambda = calclambdafromdist(ndist,L,prctile)
% Calculate mean lambda values from layer number distribution ndist, which 
% covers a distance L. Most likely values (mode of distribution) as well as 
% the associated uncertainties (percentiles of distribution) are calculated. 

%% Most likely mean layer thickness in batch:
% Use the mode of layer number distribution: 
[~,imax] = max(ndist(:,2));
nLayerML = ndist(imax,1);
lambda(:,1) = L/nLayerML;

%% And uncertainties:
% Percentiles of probability distribution:
nLayerPrctile = prctileofprobdist(ndist(:,1),ndist(:,2),prctile);

% Corresponding mean layer thicknesses:
lambda(2:1+length(prctile)) = L./nLayerPrctile;
end