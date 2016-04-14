function out = comparefloats(x,y,eps)

%% out = comparefloats(x,y,eps)
% Compare the equality of the two numbers, x and y, which may have floating 
% values (i.e. doubles/singles). The two numbers should be equal within the
% allowed rounding error "eps". 
% Copyright (C) 2015  Mai Winstrup

%% Set default value of allowed error:
if nargin<3; eps = 10^-6; end

%% Test for equality: 
if abs(x-y) < eps
    out = true;
else
    out = false;
end