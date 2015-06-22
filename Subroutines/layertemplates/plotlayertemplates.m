function hfig = plotlayertemplates(Template,TemplateInfo,Model,hfig,color,filename)

%% hfig = plotlayertemplates(Template,TemplateInfo,Model,hfig,color,filename)
% Plot layer templates, and save figure under filename. Color should be a
% (Model.order+1)x3 array with color entries corresponding to mean signal 
% (line 1) and trajectores (subsequent lines). Slightly darker colors are
% used to plot the polynomial approximations to the templates. 

% Copyright (C) 2015  Mai Winstrup

%% Figure initialization: 
% Only done if figure doesn't exist from previously.
if ~isgraphics(hfig)
    hfig = figure; 
    set(hfig,'paperunits','centimeters','papertype','A4','paperorientation',...
        'portrait','paperposition',[1 1 Model.nSpecies*7 10]);
else
    % Plot in existing figure:
    figure(hfig)
end

%% Length of 1 year:
if ~isempty(TemplateInfo)
    L = length(TemplateInfo(1).meansignal);
else
    L = 64;
end
x = 1/(2*L):1/L:1;

%% Plot mean signal and principal components:
for j = 1:Model.nSpecies
    %% Mean signal:
    subplot(2,Model.nSpecies,j); 
    if ~isempty(TemplateInfo)
        plot(x,TemplateInfo(j).meansignal,'-','color',color(1,:),'linewidth',1)
        hold on
    end
    % Polynomial approximation:
    plot(x,polyval(Template(j).mean,x),'-','color',color(1,:)*0.8,'linewidth',1.5)
    hold on
    % Add title:
    if strcmp(Model.icecore,'SyntheticData')
        title(['Species #' Model.species{j}],'fontweight','bold','interpreter','none')    
    else
        title(Model.species{j},'fontweight','bold','interpreter','none')
    end
    if j == 1; ylabel('Mean signal','fontweight','bold'); end
    
    %% Principal components:
    subplot(2,Model.nSpecies,Model.nSpecies+j);
    if ~isempty(TemplateInfo)
        for i = 1:Model.order
            plot(x,TemplateInfo(j).pc(:,1:Model.order),'-','color',color(i+1,:),'linewidth',1)
            hold on
        end
    end
    % Polynomial approximation:
    for i = 1:Model.order
        plot(x,polyval(Template(j).traj(:,i),x),'-','color',color(i+1,:)*0.8,'linewidth',1.5)
        hold on
    end
    % Make legend:
    for i = 1:Model.order
        legendname{i,:} = ['PC' num2str(i)];
    end
    legend(legendname)
    legend('boxoff')
    if j == 1; ylabel('Principal Components','fontweight','bold'); end
end

%% Save plot:
if ~isempty(filename)
    print(filename,'-djpeg','-r400')
end