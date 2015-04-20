function Z = designmatrix(Model,Template,d)

%% DESIGNMATRIX(Model,Layershape,d):
% Constructing the design matrix describing an annual layer of duration d. 
% Parts of the design matrix which correspond to unused data series are
% removed. Parameters are tied between observations and their derivatives.
% Derivatives are measured per pixel.
% Mai Winstrup 2014-03-31 17:12

% Copyright (C) 2015  Mai Winstrup

%% Constructing symmetrical layers;
x = 1/(2*d):1/d:1;

nDeriv = Model.derivatives.nDeriv;
Z = zeros((nDeriv+1)*d,Model.order,Model.nSpecies);

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



%if (strcmp(Model.type, 'PCA'))
%         Z = zeros(3*d,Model.order,Model.nSpecies);
%         for j = 1:Model.nSpecies
%             for i = 1:Model.order
%                 Z(1:d,i,j) = polyval(Template(j).traj(i,:),x);
%                 Z(d+1:2*d,i,j) = polyval(Template(j).dtraj(i,:),x)/d;
%                 Z(2*d+1:3*d,i,j) = polyval(Template(j).d2traj(i,:),x)/d^2;
%              end
%         end
    
%elseif (strcmp(Model.type, 'FFTcomp'))
%        Z = zeros(3*d,Model.order,Model.nSpecies);
        % DC component:
%        Z(1:d,1,1:Model.nSpecies) = 1;
%        Z(d+1:2*d,1,1:Model.nSpecies) = 0;
%        Z(2*d+1:3*d,1,1:Model.nSpecies) = 0;
            
%        for i = 1:(Model.order-1)/2
            % Cosine components:
%            Z(1:d,(i-1)*2+2,1:Model.nSpecies) = repmat(cos(2*pi*i*x),[1 1 2]);
%            Z(d+1:2*d,(i-1)*2+2,1:Model.nSpecies) = repmat(-(2*pi*i)*1/d*sin(2*pi*i*x),[1 1 2]);
%            Z(2*d+1:3*d,(i-1)*2+2,1:Model.nSpecies) = repmat(-(2*pi*i)^2*1/d^2*cos(2*pi*i*x),[1 1 2]);
                
            % Sine components:
 %           Z(1:d,(i-1)*2+3,1:Model.nSpecies) = repmat(sin(2*pi*i*x),[1 1 2]);
 %           Z(d+1:2*d,(i-1)*2+3,1:Model.nSpecies) = repmat((2*pi*i)*1/d*cos(2*pi*i*x),[1 1 2]);
 %           Z(2*d+1:3*d,(i-1)*2+3,1:Model.nSpecies) = repmat(-(2*pi*i)*1/d^2*sin(2*pi*i*x),[1 1 2]);
 %       end
%end

%% Removing entries corresponding to unused data series:
% if d>0
%     % Removing unused data series:
%     Zmask = true(d,nDeriv);
%     Z = Z(Zmask(:),:,:);
% end