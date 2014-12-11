function [xround ndigit]= roundsignificant(x,c)

%% ROUNDSIGNIFICANT(x,c)
% Rounding the value of x to c significant digits. Also the actual number
% of digits of the result is given (ndigit).
% Mai Winstrup, 2011

%% Rounding:
n=floor(log10(x));
% Number of actual digits in xround:
ndigit = nan(length(x),1);
for i = 1:length(x)
    ndigit(i) = max(0,-n(i)+c-1);
end

xround = round((x*10^(c-1))./10.^n).*10.^n/(10^(c-1));
% Making sure it is perfectly rounded:
xround = round(xround.*10.^ndigit)./10.^ndigit;