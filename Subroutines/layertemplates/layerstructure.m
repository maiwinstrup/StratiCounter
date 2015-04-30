function [Template, TemplateInfo] = ...
    layerstructure(data,depth,layercounts,unc_layercounts,Model,Runtype)

%% [Template, TemplateInfo] = ...
%    layerstructure(data,depth,layercounts,unc_layercounts,Model,Runtype)
% Calculating a template for an annual layer in the data file. Presently, 
% layers can only be characterized by their principal components, other 
% options may be added later. (A template based on FFT components is partly 
% implemented). Layers containing NaNs are removed from consideration when 
% calculating the principal components. 
% All layers are first linearly upsampled to contain the same number of 
% data points. All provided layer counts are used. Uncertain layers (if 
% such exist) are weighted half. 

% Copyright (C) 2015  Mai Winstrup

%% Initialize:
switch Model.type
    case 'PCA'
        Template(1:Model.nSpecies) = struct('mean',[],'dmean',[],...
            'traj',[],'dtraj',[]);
        TemplateInfo(1:Model.nSpecies) = struct('meansignal',[],...
            'pc',[],'score',[],'pcvar',[],'explained',[],'tsquare',[]); 
    case 'FFT'
        Template(1:Model.nSpecies) = struct('dc',[],'phase',[],'amplitude',[]);
        TemplateInfo(1:Model.nSpecies) = struct('meansignal',[],'meansignal_alt',[]);
end

%% Layer templates:
for j = 1:Model.nSpecies
    
    %% Stack data as layers of equal length:
    % Only original data series, not derivatives, is used for calculating
    % the layer templates. 
    stack = stacklayers(depth,data(:,1,j),layercounts,...
        unc_layercounts,Model.dtstack); 

    % Plot stack as boxplot:
    if Runtype.plotlevel > 2
        plotstack(stack,Model)
        if strcmp(Model.icecore,'SyntheticData')
            title(['Stack: Species #' Model.species{j}],'fontweight','bold')    
        else
            title(['Stack: ' Model.species{j}],'fontweight','bold')
        end
    end
    
    %% Mean annual signal:
    meansignal = nanmean(stack,1);
    
    %% Normalize data on a yearly basis before analysis?
    switch Model.normalizelayer
        case 'minmax'
            % Normalize each layer using min-max:
            minvalues = repmat(nanmin(stack,[],2),1,size(stack,2));
            maxvalues = repmat(nanmax(stack,[],2),1,size(stack,2));
            stack = (stack-minvalues)./(maxvalues-minvalues);
        
        case 'zscore'
            % Normalize each layer using zscore:
            meanvalues = repmat(nanmean(stack,2),1,size(stack,2));
            stdvalues = repmat(nanstd(stack,[],2),1,size(stack,2));
            stack = (stack-meanvalues)./stdvalues;
            
        case 'minusmean'
            % Subtract the mean signal from all layers:
            stack = stack-repmat(meansignal,size(stack,1),1);
    end
    
    % Plot normalized stack as boxplot:
    if ~strcmp(Model.normalizelayer,'none') && Runtype.plotlevel>2
        plotstack(stack,Model)
        if strcmp(Model.icecore,'SyntheticData')
            title(['Stack (' Model.normalizelayer '): Species #' num2str(j)],'fontweight','bold')    
        else
            title(['Stack (' Model.normalizelayer '): ' Model.species{j}],'fontweight','bold')
        end
    end

    %% Compute annual layer templates using one of several methods:
    switch Model.type
        case 'PCA'
            %% Principal component analysis:
            % Do not center the data first (i.e. subtract mean signal); if 
            % we wish to do so, this has been done above.
            [pc,score,pcvar,tsquare,explained] = pca(stack,'Centered',false,'NumComponents',5);
            
            % If only nans in data series:
            if isnan(pc)
                pc = nan(1,length(meansignal));
            end
            
            if Runtype.plotlevel > 1
                plotexplainedvar(explained)
            end
            
            % Keeping information corresponding to the first five
            % principal components:
            TemplateInfo(j).meansignal(:,1) = meansignal;
            TemplateInfo(j).pc = pc;
            TemplateInfo(j).score = score;
            TemplateInfo(j).pcvar = pcvar(1:min(5,length(pcvar)));
            TemplateInfo(j).explained = explained(1:min(5,length(pcvar)));
            TemplateInfo(j).tsquare = tsquare;
            
            %% Polynomial approximations of mean layer signal and the 
            % principal components:
            Template(j)=polyapprox(TemplateInfo(j).meansignal,TemplateInfo(j).pc,Model); 
            
        case 'FFT'
            %% Fourier components of layers:
            fftres = nan(size(stack,1),1/Model.dtstack);
            for i = 1:size(stack,1)
                fftres(i,:) = fft(stack(i,:));
            end
            
            % Mean DC signal:   
            a0=Model.dtstack*nanmean(fftres(:,1));
            % Mean amplitudes:
            K = 1/Model.dtstack;
            cn = fftres(:,2:K/2);
            an = nanmean(-2/K*imag(cn));
            bn = nanmean(2/K*real(cn));
            amplitude = (an.^2+bn.^2).^0.5;
    
            % Phase:
            phase = atan2(an,bn); % expansion in cosine waves
         
            % Save as template:
            Template(j).dc = a0;
            KK = (Model.order-1)/2;
            Template(j).phase = phase(1:KK);
            Template(j).amplitude = amplitude(1:KK);
            
            % Alternative mean layer shape:
            t = (1/2*K):(1/K):(1-1/K);
            meansignal_alt = nan(1,K);
            for k = 1:K
                meansignal_alt(k) = a0 + sum(amplitude(1:KK).*cos((1:KK).*2*pi*t(k)-phase(1:KK)));
            end
            TemplateInfo(j).meansignal = meansignal; % Original
            TemplateInfo(j).meansignal_alt = meansignal_alt; % Alternative based on frequency content
    end
end
end

%% Plotting subfunctions:
function plotstack(stack,Model)
figure;
boxplot(stack); 
title('Stack','fontweight','bold')
xtickloc = 0:1/(4*Model.dtstack):1/Model.dtstack;
xtickloc(1)=1;
set(gca,'xtick',xtickloc,'xticklabel',0:0.25:1)
xlabel('1 year')
end

function plotexplainedvar(explained)
figure;
plot(explained, '.-k')
xlabel('#PC')
ylabel('Explained variance')
xmax = find(explained<explained(1)/100,1,'first');
xlim([0 xmax])
title('Principal Component Analysis','fontweight','bold')
end