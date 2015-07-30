function depth = makedepthscale(dstart0,dend,dx,dx_offset)

%% depth = makedepthscale(dstart,dend,dx,dx_offset)
% Copyright (C) 2015  Mai Winstrup

%% Actual starting value should be:
dstart = floor(dstart0)+dx_offset*dx;

%% Corresponding depth scale:
depth = (dstart:dx:dend)';

% Remove parts from outside interval:
depth = depth(depth>=dstart0);