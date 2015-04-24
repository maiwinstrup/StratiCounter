function Template = polyapprox(meansignal,pc,Model)

%% Template = polyapprox(meansignal,pc,Model)
% Calculate polynomial approximations (and derivatives) corresponding to 
% the "raw" input as given by mean signal and principal components (pc).
% Copyright (C) 2015  Mai Winstrup

%% Coefficients for fitted polynomials for mean signal: 
% Mean signal:
x = (1/(2*length(meansignal)):1/length(meansignal):1);
Template.mean = polyfit(x(:),meansignal(:),Model.pcPolOrder);

% And its derivatives:
coeff = Template.mean;
for k = 1:Model.derivatives.nDeriv
    Template.dmean(:,k) = polyder(coeff);
    % Coefficients of current derivative:
    coeff = Template.dmean(:,k);
end
    
%% Similarly for the principal components: 
% Ensure format of pc matrix to be [length(x),Model.order]:
if size(pc,2)==length(x); pc = pc'; end

% Polynomial approximations:
for i = 1:Model.order
    % Signal:
    Template.traj(:,i) = polyfit(x(:),pc(:,i),Model.pcPolOrder);
    % Derivatives:
    coeff = Template.traj(:,i);
    for k = 1:Model.derivatives.nDeriv
        Template.dtraj(:,i,k) = polyder(coeff);
        coeff = Template.dtraj(:,i,k);
    end
end