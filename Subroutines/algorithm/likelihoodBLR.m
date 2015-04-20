 function [logb, ExpVal] = likelihoodBLR(data,Template,Layerpar,d,i0,T,Model)

%% [logb, E] = likelihoodBLR(data,d,Model,Template,Layerpar,i0,T):
% Calculating the log-likelihood that the observation sequence from t1 to
% t2 spans exactly one year:
% logb(t1,t2) = log(P(obs(t1:t2) | S([t1:t2])=j))
% This is calculated for all values of t1 and t1+min(d)<=t2=<t1+max(d),
% while assuming the layers to be described by an "layershape" with
% parameters layerpar.par(:,j), which have a covariance matrix layerpar.cov. 
% To this signal is added an i.i.d. gaussian noise distribution with
% variance layerpar.nvar(1:3,j). 
% Probabilities are calculated using Bayesian Linear Regression.

% Variables:
% data(:,1:3,j): array containing the observation series (and derivatives)
%   for species j
% d: possible durations
% model: e.g. type and order of model 
% layershape: shape of annual layer curve
% layerpar.par(:,j): expectancy value of parameters for species j
% layerpar.cov: covariance matrix of parameters (one matrix for all species)
% layerpar.nvar(1:3,j): spread of white noise in data series j compared to model

% Copyright (C) 2015  Mai Winstrup

%% Duration parameters:
dmin = min(d); % Minimum duration
dmax = max(d); % Maximum duration
D = length(d); % Number of possible durations

%% Appending the observation sequence(s) with NaNs in both ends, and remove 
% unused derivative data series:
obs = [NaN(dmax,Model.derivatives.nDeriv+1,Model.nSpecies); data(:,1:Model.derivatives.nDeriv+1,:); ...
    NaN(dmax,Model.derivatives.nDeriv+1,Model.nSpecies)];

%% Calculating the log-likelihood:
% For t<1, is calculated the conditional probability: 
% P(obs(1:(t+d) | S([t:t+d])=j)
% And for tend>T is calculated: P(obs(t:T) | S([t:t+d])=j)

% Inverse of parameter covariance matrix ("precision matrix"):
Pmatrix = inv(Layerpar.cov);

% Initialization of matrices: 
% Matrices are indexed by t (ending time), d (duration), and j (species).
logb = ones(T+dmax,D,Model.nSpecies)*-inf;
% E contains various expectation values for segments:
ExpVal = repmat(struct('invw',[],'N',0,'r',...
    nan(Model.order,1),'Cr',[],'rr',nan(Model.order,Model.order),...
    'X',[],'XWX',nan(Model.order,Model.order),'XWoXr',nan(Model.order,1),...
    'oXr',[]), [T+dmax, D, Model.nSpecies]);
%the coder is dumb :)
coder.varsize('ExpVal(:).invw');
coder.varsize('ExpVal(:).Cr');
coder.varsize('ExpVal(:).X');
coder.varsize('ExpVal(:).oXr');

for istart = i0(1,1):T+dmax
    for id = max(1,dmax-dmin-istart+3):D
        iend = istart+d(id)-1;
        obs_segment = obs(istart:iend,:,:);
        tend = iend-dmax;
        if sum(isfinite(obs_segment(:))) > 0
            [logb(tend,id,:), ExpVal(tend,id,:)] = ...
                calc_b(obs_segment,Model,Template,Layerpar,d(id),Pmatrix);
        else
            logb(tend,id,:) = 0;  % Corresponding to a probability of 1
        end
    end
end
end
 
function [logb, ExpVal] = calc_b(obs_segment,Model,Template,Layerpar,d,P0matrix)
%% [logb, ExpVal] = calc_b(obs_segment,model,layershape,layerpar,d,P0matrix)
% Calculating log(b) using Bayesian Linear Regression to a given template. 

%% x-values:
x = 1/(2*d):1/d:1;

%% If minusmean: Subtract mean signal from observation segment
if strcmp(Model.normalizelayer,'minusmean')
    meansignal = nan(d,Model.derivatives.nDeriv,Model.nSpecies);
    for j = 1:Model.nSpecies
        meansignal(:,1,j) = polyval(Template(j).mean,x)';
        for k = 1:Model.derivatives.nDeriv
            meansignal(:,k+1,j) = polyval(Template(j).dmean(k,:,:),x)'/d^k;
        end
        % Derivatives of mean signal are calculated per pixel distance 
        % (not per time). 
    end
    obs_segment = obs_segment - meansignal(:,1:Model.derivatives.nDeriv+1,:);
end

%% Reshape array, as we consider the signal and its derivatives to be 
% modelled with the same parameter values:
obs_segment = reshape(obs_segment,(Model.derivatives.nDeriv+1)*d,Model.nSpecies);
% Indexed: obs_segment(1:Nderiv*d,1:model.nSpecies)

%% Construct design matrix for the linear regression:
X0 = designmatrix(Model,Template,d); % 2014-03-31 17:17
% Indexed: X0(1:Nderiv*d,1:Model.order,1:Model.nSpecies)

%% Initialize matrices:
logb = nan(1,Model.nSpecies);
ExpVal = repmat(struct('invw',[],'N',0,'r',...
    nan(Model.order,1),'Cr',[],'rr',nan(Model.order,Model.order),...
    'X',[],'XWX',nan(Model.order,Model.order),'XWoXr',nan(Model.order,1),...
    'oXr',[]), [1, Model.nSpecies]);
%the coder is dumb :)
coder.varsize('ExpVal(:).invw');
coder.varsize('ExpVal(:).Cr');
coder.varsize('ExpVal(:).X');
coder.varsize('ExpVal(:).oXr');

%% Calculating logb and E for a given segment:
for j = 1:Model.nSpecies
    % Observation segment:
    mask = isfinite(obs_segment(:,j));
    obs_j = obs_segment(mask,j);
    % Design matrix for the linear regression for species j, where data is
    % available:
    X = X0(mask,:,j);

    % Form diagonal invW matrix, controlling the weighting of white noise
    % variance component on derivative data series:
    w = ones(d,1)*Model.derivnoise;
    w = w(:); 
    w = w(mask);
    W = diag(w);
    % And its inverse:
    ExpVal(:,j).invw = 1./w;
    invW = diag(ExpVal(:,j).invw);
    
    % Parts of X and W corresponding to unused data series and nans in data
    % have been removed.
    % Total number of used data points:
    ExpVal(:,j).N = sum(mask);
    
    % Selecting the part of cov and P matrices corresponding to this
    % species. Individual species are considered uncorrelated. 
    PHI = Layerpar.cov((j-1)*Model.order+1:j*Model.order,(j-1)*Model.order+1:j*Model.order);
    P = P0matrix((j-1)*Model.order+1:j*Model.order,(j-1)*Model.order+1:j*Model.order);
    
    %% Bayesian Linear Regression: 
    obsj_mean = X*Layerpar.par(:,j);
    Shat = X*PHI*X'+Layerpar.nvar(j)*W; 
    
    % Resulting probability density of data segment:
    L = chol(Shat,'lower');
    v = L\(obs_j-obsj_mean);
    logb(j) = -ExpVal(:,j).N/2*log(2*pi)-sum(log(diag(L)))-0.5*sum(v.^2);
    
    %% Expectation values etc:
    % Conditioned on the observations and the present set of parameters.
    % Expectation value of rj, conditioned on oj:
    M = X'*invW*X+Layerpar.nvar(j)*P;
    ExpVal(:,j).r = M\(X'*invW*(obs_j-obsj_mean));
    % Covariance of rj, conditioned on oj: 
    ExpVal(:,j).Cr = Layerpar.nvar(j)*eye(Model.order,Model.order)/M;
    % Expectation value of rr', conditioned on oj (for computing "cov"):
    ExpVal(:,j).rr = ExpVal(:,j).r*ExpVal(:,j).r'+ ExpVal(:,j).Cr;

    % Matrices for computing "par":     
    ExpVal(:,j).X = X;
    ExpVal(:,j).XWX = X'*invW*X;
    ExpVal(:,j).XWoXr = X'*invW*(obs_j-X*ExpVal(:,j).r);

    % Matrices for computing "nvar":
    ExpVal(:,j).oXr = obs_j-X*ExpVal(:,j).r; % An estimate of X*beta
end
end