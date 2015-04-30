function Z = designmatrix(Template,d,Model)

%% Z = designmatrix(Template,d,Model)
% Constructing the design matrix describing an annual layer of duration d. 
% Parameters are tied between observations and their derivatives. 
% Derivatives are measured per pixel.
% Copyright (C) 2015  Mai Winstrup

%% Constructing symmetrical layers;
x = 1/(2*d):1/d:1;

nDeriv = Model.derivatives.nDeriv;
Z = zeros((nDeriv+1)*d,Model.order,Model.nSpecies);

if strcmp(Model.type, 'PCA')
    % Trajectories:
    for j = 1:Model.nSpecies
        for i = 1:Model.order
            Z(1:d,i,j) = polyval(Template(j).traj(:,i),x);
        end
    end

    % Derivatives:
    for j = 1:Model.nSpecies
        for i = 1:Model.order
            for k = 1:nDeriv
                Z(k*d+1:(k+1)*d,i,j) = polyval(Template(j).dtraj(:,i,k),x)/d^k;
            end
        end
    end
    
elseif strcmp(Model.type, 'FFT')
    % DC component:
    Z(1:d,1,:) = 1;
    % Derivatives are equal to 0:
    Z(d+1:(nDeriv+1)*d,1,:) = 0;
    
    for i = 1:(Model.order-1)/2
        % Cosine component:
        Z(1:d,(i-1)*2+2,:) = repmat(cos(2*pi*i*x),[1 1 nDeriv]);
        % Sine component:
        Z(1:d,(i-1)*2+3,:) = repmat(sin(2*pi*i*x),[1 1 nDeriv]);
        
        for k = 1:2:nDeriv
            % Cosine component:
            Z(k*d+1:(k+1)*d,(i-1)*2+2,:) = ...
                repmat(-(2*pi*i)^k*1/d^k*sin(2*pi*i*x),[1 1 nDeriv]);
            % Sine components:
            Z(k*d+1:(k+1)*d,(i-1)*2+3,:) = ...
                repmat(-(2*pi*i)^k*1/d^k*cos(2*pi*i*x),[1 1 nDeriv]);
        end
        
        for k = 2:2:nDeriv
            % Cosine component:
            Z(k*d+1:(k+1)*d,(i-1)*2+2,:) = ...
                repmat(-(2*pi*i)^k*1/d^k*cos(2*pi*i*x),[1 1 nDeriv]);
            % Sine components:
            Z(k*d+1:(k+1)*d,(i-1)*2+3,:) = ...
                repmat(-(2*pi*i)^k*1/d^k*sin(2*pi*i*x),[1 1 nDeriv]);
        end
    end
end