function data_new = boxcoxtransform(data,lambda,alpha,plotlevel)

%% data_new = boxcoxtransform(data,lambda,alpha,plotlevel)
% Performing a box-cox transformation on data, with power parameter lambda
% and shift parameter alpha. (Strictly speaking: In a box-cox
% transformation, the geometric mean is assumed equal to 1).

% Copyright (C) 2015  Mai Winstrup
% 2014-05-18 18:47: Minor adjustments
% 2014-07-16 09:24: Showplots -> plotlevel

%% Box-Cox transformation:
% Geometric mean of intensity values:
gm = geomean(data(isfinite(data)));
if lambda == 0
    data_new = gm*log(data+alpha);
elseif lambda == 1
    data_new = data+alpha-1;
else
    data_new = ((data+alpha).^lambda-1)/(lambda*gm^(lambda-1));
end
   
%% Plot transformed data:
if plotlevel > 0
    figure;
    subplot(2,1,1)
    plot(data,'-k')
    title('Original data','fontweight','bold')
    subplot(2,1,2)
    plot(data_new,'-b')
    title(['Box-Cox transformated data (\lambda=' num2str(lambda) ', \alpha=' num2str(alpha) ')'],'fontweight','bold')
end