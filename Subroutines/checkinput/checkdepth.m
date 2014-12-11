function newinterval = checkdepth(interval,dx,dx_center)

%% newinterval = checkdepth(interval,dx,dx_center)
% This function checks that the boundaries of the given interval do indeed 
% coincide with the depth of a specific datapoint. If not, newinterval is 
% the slightly smaller interval for which this is the case. The depth scale 
% is determined based on dx and dx_center.

% Mai Winstrup,
% 2014-08-22 15:14

%% Calculate new depths:
newinterval(1) = ceil((interval(1)-dx_center)/dx)*dx+dx_center;
newinterval(2) = floor((interval(2)-dx_center)/dx)*dx+dx_center;