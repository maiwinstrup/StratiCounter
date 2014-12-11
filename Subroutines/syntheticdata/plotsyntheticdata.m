function plotsyntheticdata(data,depth,target,par,dur_px,counts,SynthData,filename,plotlevel)

%% plotsyntheticdata(data,depth,target,par,dur_px,counts,SynthData,filename,plotlevel)
% Plot target data and noisy observations, and distributions of layer 
% thicknesses and model parameters (input and output values).

% Mai Winstrup
% 2014-06-20 23:08: Separate script
% 2014-08-19 20:56: Changes to when to plot, and save figure

%% Target data and noisy observations:
figure;
for j = 1:SynthData.nSpecies
    subplot(SynthData.nSpecies,1,j)
    hline(1)=plot(depth,target(:,j), '-k','linewidth',1.1);
    hold on
    hline(2)=plot(depth,data(:,1,j), '-b','linewidth', 1.2); % With added noise and then smoothed 
    yvaluelim = quantile(data(:,1,j),[0.05 0.95]);
    for i = 1:30
        plot(counts(i)*[1 1], yvaluelim, '-r')
    end
    ylabel(['Data series #' num2str(j)])
    xlim([0 25*mean(diff(counts))])
    if j == 1
        hleg = legend(hline,{'Target','Observed'});
        set(hleg,'location','northwest')
        legend boxoff
    end
end
suptitle('Target signal and observations')
xlabel('Depth [m]')

% Save figure:
if ~isempty(filename)
    print(filename,'-djpeg','-r400')    
end
% And maybe close it:
if plotlevel==1; close(gcf); end

%% Layer thickness distribution:
% Observed layer thicknesses and corresponding values of my and sigma:
lambda_out = dur_px*SynthData.dx;
my_out = nanmean(log(lambda_out));
sigma_out = nanstd(log(lambda_out));
   
disp('Input values of my and sigma:')
disp(['my = ' num2str(SynthData.Modelpar.my)]);
disp(['sigma = ' num2str(SynthData.Modelpar.sigma)]);
  
disp('Values of my and sigma from actual distribution:')
disp(['my = '  num2str(my_out)])
disp(['sigma = ' num2str(sigma_out)])

%% Plot input and output distributions: 
if plotlevel>1
    % Layer thickness distribution:
    figure;
    dval = min(lambda_out)-SynthData.dx:SynthData.dx:max(lambda_out)+SynthData.dx;
    hist(lambda_out,dval)
    hold on
    nLayers = length(counts);
    plot(dval,lognpdf(dval,my_out,sigma_out)*nLayers*SynthData.dx,'-g','linewidth',2)
    plot(dval,lognpdf(dval,SynthData.Modelpar.my,SynthData.Modelpar.sigma)*nLayers*SynthData.dx,'--r','linewidth',2)
    legend({'Actual distribution','Log-normal fit to actual distribution','Log-normal input distribution'})
    title('Layer thickness distribution','fontweight','bold')
        
    %% Layer parameter value distributions: 
    cov_values = reshape(diag(SynthData.Modelpar.cov),SynthData.order,SynthData.nSpecies);
    mean_par_out = nanmean(par);
    var_par_out = nanvar(par);
    
    for j = 1:SynthData.nSpecies
        xmin = SynthData.Modelpar.par(:,j)-3*cov_values(:,j).^0.5;
        xmax = SynthData.Modelpar.par(:,j)+3*cov_values(:,j).^0.5;
        dx = (xmax-xmin)/30;

        figure;
        for i = 1:SynthData.order
            subplot(1,SynthData.order,i)
            x = xmin(i):dx(i):xmax(i);
            hist(par(:,i,j),x)
            title(['par(' num2str(i) ')'],'fontweight','bold')
            % Fitted distribution
            dist = normpdf(x,mean_par_out(1,i,j),var_par_out(1,i,j)^0.5);
            hold on
            A = nLayers*dx(i);
            plot(x,dist*A,'-g','linewidth',2)
            % And original distribution:
            dist0 = normpdf(x,SynthData.Modelpar.par(i,j),cov_values(i,j)^0.5);
            plot(x,dist0*A,'--r','linewidth',2)           
        end
        suptitle(['Layer parameters, data record #' num2str(j)])
    end
end