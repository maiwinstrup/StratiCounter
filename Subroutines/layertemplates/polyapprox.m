function Template = polyapprox(meansignal,pc,Model)

%TODO: Laura needs to write a C-version of this function, as it's being a
%jerk about conversion.

%% Template = polyapprox(meansignal,pc,Model)
% Calculate polynomial approximations (and derivatives) corresponding to 
% the "raw" input as given by mean signal and principal components (pc).

% Copyright (C) 2015  Mai Winstrup
% 2014-08-17 17:51: Initial version
% 2014-08-20 16:28: x and meansignal -> x(:), meansignal(:)
% 2014-08-21 14:08: Ensure correct format of pc

%% Coefficients for fitted polynomials for mean signal: 
x = (1/(2*length(meansignal)):1/length(meansignal):1);
Template.mean = polyfit(x(:),meansignal(:),Model.pcPolOrder);
% And its derivatives:
Template.dmean = polyder(Template.mean);
Template.d2mean = polyder(Template.dmean);
    
%% Similarly for the principal components: 
% Ensure correct format of pc matrix (namely [length(x),Model.order]):
if size(pc,2)==length(x); pc = pc'; end

% Polynomial approximations:
for i = 1:Model.order
    Template.traj(i,:) = polyfit(x(:),pc(:,i),Model.pcPolOrder);
    Template.dtraj(i,:) = polyder(Template.traj(i,:));
    Template.d2traj(i,:) = polyder(Template.dtraj(i,:));
end