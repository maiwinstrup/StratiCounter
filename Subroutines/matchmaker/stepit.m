function outdata = stepit(indata)
% STEPIT.m : produces step curve data from ordinary [depth data1 data2 ...] array. 
% Asumes that the depth given for each data value is the end depth
% (that is : depth of deepest end) of the section to which the measured
% value correspond.
% Copyright (C) Sune Olander Rasmussen

[L, N] = size(indata);

if N < 2
    disp('ERROR in stepit.m : data must have at least two columns');
    outdata = [];
    return;
end;

outdata = zeros(2*L, N);
outdata(:,1) = indata([ 1 ceil(0.5 : 0.5 : L-0.5)],1);
outdata(1,1) = indata(1,1) - (indata(2,1) - indata(1,1));
outdata(:,2:N) = indata(ceil(0.5 : 0.5 : L),2:N);