function plotlayercounts(counts,data)
% Copyright (C) 2015  Mai Winstrup

% Plot layer counts onto figure:
hbar0 = quantile(data,0.05);
hbar1 = quantile(data,0.95);
hbar2 = quantile(data,0.90);  
if ~isempty(counts)
    mask = counts(:,3)==0;
    plot(counts(mask,1)*[1 1],[hbar0 hbar1],'-r')
    if sum(counts(:,3))>0
       plot(counts(~mask,1)*[1 1],[hbar0 hbar2],'--r')
    end
end