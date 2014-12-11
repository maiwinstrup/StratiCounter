function [xstep, ystep]=stepit(x,y)
% A function to step data for plotting purposes

dx = (diff(x))';
x = [x(1)-mean(dx)/2,x(2:end)'-dx/2; x(1:end-1)'+dx/2, x(end)+mean(dx)/2];
xstep = x(:);

y = [y(:)'; y(:)'];
ystep = y(:);

