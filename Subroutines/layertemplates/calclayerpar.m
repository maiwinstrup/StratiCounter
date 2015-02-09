function [layerpar, ParML, ParMAP, Layerpar0] = calclayerpar(Model,Data,layerpos,unc,Template,Runtype)

%% [layerpar,ParML, ParMAP, Layerpar0] = calclayerpar(Model,Data,layerpos,unc,Template,Runtype)
% Computing maximum likelihood Model parameters for each layer individually, 
% based on manually counted layer positions. 
% If some layer boundaries are given as uncertain, we have:
% 1st column: Composite of the below two data sets
% 2nd column: All layers are considered as certain
% 3rd column: All uncertain layers are removed
% An easy way to weight the layer parameters according to the certainty of 
% layer boundaries is to use the composite of both data sets (first row).
% Based on the individual layerparameter values, also the maximum 
% likelihood estimate of the model parameters are calculated. 

% Copyright (C) 2015  Mai Winstrup
% 2014-10-01 16:43: New version

%% Layer depths: 
% In case of uncertain layer boundaries in data set, the array contains
% both two separate data series (with/without uncertain layers). This is
% done to ensure that also inference regarding the relationship between
% successive layers can be made. Furthermore, a composite data set is
% configured. If no uncertain layer boundaries, only one set of layer
% parameter values are made. 
if sum(unc)>0
    layerpos1 = layerpos;
    layerpos2 = layerpos;
    layerpos2(unc==1)=[];
    layerpos = [layerpos1; NaN; layerpos2];
end

% Depth of start of layers:
if sum(unc)==0
    layerpar(1).depth = layerpos(1:end-1);
else
    layerpar(2).depth = layerpos1(1:end-1);
    layerpar(3).depth = layerpos2(1:end-1);
    [layerpar(1).depth, index] = ...
        sort([layerpar(2).depth; layerpar(3).depth],1,'ascend');
    % The remaining layer parameters are sorted similarly.
end

%% Annual layer thicknesses:
% If the data series has uncertain manually counted layer boundaries, we 
% use two sets of data: One in which all layer boundaries are assumed
% certain, and one where all uncertain layer boundaries have been removed.
% The composite of these two data sets is then analyzed to obtain ParML.
lambda = layerpos(2:end)-layerpos(1:end-1); %[m]
lambda = lambda(isfinite(lambda));
    
if sum(unc)==0
    layerpar(1).lambda = lambda;
else
    N = length(layerpos1)-1;
    layerpar(2).lambda = lambda(1:N);
    layerpar(3).lambda = lambda(N+1:end);
    lambda_tot = [layerpar(2).lambda; layerpar(3).lambda];
    layerpar(1).lambda = lambda_tot(index);
end

%% Converting layer boundaries to pixels: 
T = size(Data.data,1);
if ~isempty(Model.dx)
    % layerpos_px is the first pixel of layer: 
    layerpos_px = interp1([Data.depth(1)-Model.dx;Data.depth;Data.depth(end)+Model.dx],0:T+1,layerpos+Model.dx/2,'nearest',nan);
else
    dx1 = Data.depth(2)-Data.depth(1);
    dx2 = Data.depth(end)-Data.depth(end-1);
    layerpos_px = interp1([Data.depth(1)-dx1;Data.depth;Data.depth(end)+dx2],0:T+1,layerpos+Model.dx/2,'nearest',nan);
end
layerpos_px(isnan(layerpos(:)))=nan;
M = length(layerpos_px);

%% Layer parameter estimates from each individual layer:
for j = 1:Model.nSpecies
    par_hat=nan(Model.order,M-1);
    nvar_hat=nan(length(Model.deriv),M-1);
    eps = cell(1,M-1);
 
    for i=1:M-1
        if isfinite(layerpos_px(i))&&isfinite(layerpos_px(i+1))
            % Selecting appropriate data series (species j):
            datasegment=Data.data(layerpos_px(i):layerpos_px(i+1)-1,Model.deriv+1,j);
            d=size(datasegment,1);
        
            if strcmp(Model.normalizelayer,'minusmean')
                % Subtracting mean signal from data:
                x = 1/(2*d):1/d:1;
                clear datasegment_subtract
                datasegment_subtract(:,1) = polyval(Template(j).mean,x)';
                datasegment_subtract(:,2) = polyval(Template(j).dmean,x)'/d;
                datasegment_subtract(:,3) = polyval(Template(j).d2mean,x)'/d^2;
                datasegment = datasegment - datasegment_subtract(:,Model.deriv+1);
            end
            mask=isfinite(datasegment(:));
            datasegment = datasegment(mask);
 
            if sum(mask)>Model.order
                X = designmatrix(Model,Template,d);
                % Removing parts of X-matrix corresponding to other species and
                % to NaNs in current data series:
                X = X(mask,:,j);
            
                % MLE of best fitting parameter values:
                par_hat(:,i) = X\datasegment;
            
                % Corresponding residuals:
                % Divided into residuals of the data series and it derivatives.
                data_hat = X*par_hat(:,i);
                res = nan(1,d*length(Model.deriv));
                res(mask) = datasegment-data_hat;
                res = reshape(res,d,length(Model.deriv));
                eps{1,i} = res;
                for m = 1:length(Model.deriv)
                    nvar_hat(m,i) = nanmean(res(:,m).^2);
                end
               
                if Runtype.plotlevel >= 2 && i<3
                    plotlayerfit(datasegment,datasegment_subtract,X,par_hat(:,i),d,mask,j,Model)
                end
            end
        end
    end
    
    %% Add to layerpar array: 
    if sum(unc)==0
        layerpar(1).par(:,:,j) = par_hat;
    else
        layerpar(2).par(:,:,j) = par_hat(:,1:N);
        layerpar(3).par(:,:,j) = par_hat(:,N+1:end);
        par_tot = [layerpar(2).par(:,:,j), layerpar(3).par(:,:,j)];
        layerpar(1).par(:,:,j) = par_tot(:,index);
    end

    % Mean of white noise variance:
    if sum(unc)==0
        layerpar(1).nvar(:,:,j) = nvar_hat;
    else
        layerpar(2).nvar(:,:,j) = nvar_hat(:,1:N);
        layerpar(3).nvar(:,:,j) = nvar_hat(:,N+1:end);
        nvar_tot = [layerpar(2).nvar(:,:,j) layerpar(3).nvar(:,:,j)];
        layerpar(1).nvar(:,:,j) = nvar_tot(:,index);
    end

    % White noise:
    if sum(unc)==0
        layerpar(1).eps(:,j) = eps;
    else
        layerpar(2).eps(:,j) = eps(1:N);
        layerpar(3).eps(:,j) = eps(N+3:end);
        eps_tot = [layerpar(2).eps(:,j); layerpar(3).eps(:,j)];
        layerpar(1).eps(:,j) = eps_tot(index);
    end
end

%% Maximum Likelihood estimmate of annual layer parameters
% This estimate is based on the composite set.
mask = layerpar(1).depth>=Model.initialpar(1) & layerpar(1).depth<Model.initialpar(2);
ParML.my = nanmean(log(layerpar(1).lambda(mask)));
ParML.sigma = nanstd(log(layerpar(1).lambda(mask)));

ParML.cov = zeros(Model.nSpecies,Model.nSpecies);
for j = 1:Model.nSpecies
    ParML.par(:,j) = nanmean(layerpar(1).par(:,mask,j),2);
    ParML.cov(Model.order*(j-1)+1:Model.order*j,Model.order*(j-1)+1:Model.order*j) = nancov(layerpar(1).par(:,mask,j)');
    ParML.nvar(:,j) = nanmean(layerpar(1).nvar(:,mask,j),2);
end

%% Setting weight for the white noise component of derivative data series:
if strcmp(Model.wWhiteNoise,'manual')
    % Using the derived ML value for the w-values:
    Model.wWhiteNoise = nan(length(Model.deriv),Model.nSpecies);
    Model.wWhiteNoise(Model.deriv+1,:)= ParML.nvar./repmat(ParML.nvar(1,:),2,1);
    ParML.nvar = ParML.nvar(1,:);
else
    % If using analytical values:
    for j = 1:Model.nSpecies
        ParML.nvar(1,j) = 1/length(Model.deriv)*...
            sum(1./Model.wWhiteNoise(Model.deriv+1)'.*ParML.nvar(:,j));
    end
    ParML.nvar(2,:) = [];
end

%% Using the EM-algorithm to find the MAP estimate for layer parameters:
% Number of iterations:
nIter = 10;

for j = 1:Model.nSpecies
    % Initialization:
    Err = nan(Model.order,Model.order,M-1);
    XWX = nan(Model.order,Model.order,M-1);
    design = cell(M-1,1);
    XWoXr = nan(Model.order,M-1);
    oXr = cell(M-1,1);
    sumMd = nan(M-1,1);
    Par(1,nIter) = struct('par',[],'cov',[],'nvar',[]);
    invW2 = cell(M-1,1);

    % Initial parameter input:
    Par(1).par = ParML.par(j);
    Par(1).cov = ParML.cov(Model.order*(j-1)+1:Model.order*j,Model.order*(j-1)+1:Model.order*j);
    Par(1).nvar = ParML.nvar(j);

    for k = 1:nIter
        for i=1:M-1
            if isfinite(layerpos_px(i))&&isfinite(layerpos_px(i+1))
                % Picking appropriate data segment and data series:
                datasegment=Data.data(layerpos_px(i):layerpos_px(i+1)-1,1+Model.deriv,j);
                d=size(datasegment,1);
            
                if strcmp(Model.normalizelayer,'minusmean')
                    % Subtracting mean signal from data:
                    x = 1/(2*d):1/d:1;
                    clear datasegment_subtract
                    datasegment_subtract(:,1) = polyval(Template(j).mean,x)';
                    datasegment_subtract(:,2) = polyval(Template(j).dmean,x)'/d;
                    datasegment_subtract(:,3) = polyval(Template(j).d2mean,x)'/d^2;
                    datasegment = datasegment - datasegment_subtract(:,Model.deriv+1);
                end
                mask=isfinite(datasegment(:));
                datasegment = datasegment(mask);

                if sum(mask)>Model.order
                    % Design matrix:
                    X = designmatrix(Model,Template,d);
                    X = X(mask,:,j);
            
                    % Relative white noise levels:
                    diagonalvector = [Model.wWhiteNoise(1)*ones(1,d) Model.wWhiteNoise(2)*ones(1,d) Model.wWhiteNoise(3)*ones(1,d)];
                    invW = diag(diagonalvector(mask).^-1);
            
                    % Expectation value of random component r:
                    Er = (Par(k).nvar*eye(Model.order,Model.order)/Par(k).cov(:,:)+X'*invW*X)\(X'*invW*(datasegment-X*Par(k).par));
                    % Covariance of random component r:
                    Cr = Par(k).nvar*eye(Model.order,Model.order)/(Par(k).nvar*eye(Model.order,Model.order)/Par(k).cov + X'*invW*X);
                    % Conditional expectation value of r*r' for this layer:
                    Err(:,:,i)=Er*Er'+Cr;
                
                    % Value of X and X'*W^-1*X for layer:
                    design{i}=X;
                    XWX(:,:,i) = X'*invW*X;
    
                    % Value of X'*invW*(O-X*Er) for layer:
                    XWoXr(:,i) = X'*invW*(datasegment-X*Er);
                
                    % Value of the residuals O-Xr for layer:
                    oXr{i} = datasegment-X*Er;
    
                    % Number of observations:
                    sumMd(i) = length(datasegment); 
    
                    % Trace of matrix:
                    tr(i) = trace(XWX(:,:,i)*Cr);
                    invW2{i} =invW;
                end
            end
        end

        % Updated parameter value:
        Par(k+1).par = nansum(XWX,3)\nansum(XWoXr,2);  % obs- nan?
    
        % Updated covariance matrix:
        Par(k+1).cov = nanmean(Err,3);
 
        % Expectation value of E*W^-1*E:
        EE = nan(M-1,1);
        for i = 1:M-1
            if ~isempty(oXr{i})
                EE(i) = (oXr{i}-design{i}*Par(k+1).par)'*invW2{i}*(oXr{i}-design{i}*Par(k+1).par)+tr(i);
            end
        end
    
        Par(k+1).nvar = nansum(EE)/nansum(sumMd);
    end

    if Runtype.plotlevel>=2
        plotparameteriterations(Par,Model,nIter)
    end
    
    ParMAP.par(j) = Par(nIter+1).par;
    ParMAP.cov(Model.order*(j-1)+1:Model.order*j,Model.order*(j-1)+1:Model.order*j) = Par(nIter+1).cov;
    ParMAP.nvar(j) = Par(nIter+1).nvar;
end

%% Using this as initial Layerparameter estimate:
Layerpar0.my = ParML.my;
Layerpar0.sigma = ParML.sigma;
Layerpar0.par = ParMAP.par;
Layerpar0.cov = ParMAP.cov; 
Layerpar0.nvar = ParMAP.nvar; % Denne er kun en enkelt v�rdi?? JA! 

end

%% Embedded plotting functions:
function plotlayerfit(datasegment,datasegment_subtract,X,par_hat,d,mask,j,Model)   
% Plot layer data and the obtained fit to layer. 

fitresult = nan(1,d*length(Model.deriv));
fitresult(mask) = X*par_hat;
fitresult = reshape(fitresult,d,length(Model.deriv));
segment_plot = nan(1,d*length(Model.deriv));
segment_plot(mask)=datasegment;
segment_plot = reshape(segment_plot,d,length(Model.deriv));

figure;
clf
if ismember(Model.type,'PCA')
    for k = 1:length(Model.deriv)
        subplot(length(Model.deriv),1,k)
        if strcmp(Model.normalizelayer,'minusmean') 
            plot(segment_plot(:,k)+datasegment_subtract(:,Model.deriv(k)+1),'-k','linewidth',2)
            hold on
            plot(fitresult(:,k)+datasegment_subtract(:,Model.deriv(k)+1),'-b','linewidth',2)
        else
            plot(segment_plot(:,k),'-k','linewidth',2)
            hold on
            plot(fitresult(:,k),'-b','linewidth',2)
        end
        title(['Data series #' num2str(Model.deriv(k))])
    end
    suptitle(['Fit to manual layers: ' Model.species{j}])    
end
end

%% Plotting the evolution of parameters with iterations:
function plotparameteriterations(Par,Model,nIter)
% Model parameters are:
modelpar = [Par(1,:).par];
c = reshape([Par(1,:).cov],Model.order,Model.order,nIter+1);
nvar = [Par(1,:).nvar];
    
% Evolution of par:
figure;
for i = 1:Model.order
    subplot(Model.order,1,i)
    plot(modelpar(i,:))
end
xlabel('# Iteration')
suptitle('par')
    
% Evolution of indices in covariance matrix:
figure;
for i = 1:Model.order^2
    subplot(Model.order,Model.order,i)
    jj = mod(i-1,Model.order)+1;
    j2 = floor((i-1)/Model.order)+1;
    x = c(j2,jj,:);
    plot(x(:))
    if j2 == Model.order; xlabel('# Iteration'); end
end
suptitle('cov')    
    
% Evolution of nvar:
figure;
plot(nvar)
title('nvar')
xlabel('# Iteration')
end