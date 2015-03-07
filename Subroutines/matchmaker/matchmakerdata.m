function matchmakerdata(species,Model,datafilename)
%% matchmakerdata(species,Model,datafilename)
% Save unprocessed datafiles in matchmaker format. 

nSpecies = length(species);

data0 = loadrawdata(species{1},Model.path2data);
% Data:
data{1}=data0(:,2);
species{1}=species{1};
% Depth:
depth{1}=data0(:,1); 
depth_no(1)=1;
    
% Subsequent data files:    
for j = 2:nSpecies
    data0 = loadrawdata(species{j},Model.path2data);

    % Data: 
    data{j}=data0(:,2);
        
    % Depth scale:
    % On same depthscale as previous data files?
    for i = 1:length(depth)            
        if isequal(data0(:,1),depth{i})
            % Same depth scale:
            depth_no(j)=i;
            flag = 1; 
            break % Do not check further depth scales
        else
            flag = 0; 
        end
    end
    % Use old depthscale:
    if flag == 1
        continue % Continue to next data file
    end
     
    % New depth scale:
    nd = length(depth)+1; % Current "depth scale number"
    depth{nd}=data0(:,1);
    depth_no(j)=nd;
end
        
% Species are given random colours:
colours = rand(nSpecies,3);   
mask = colours>0.9;
colours(mask)=0.7*colours(mask); % avoiding too bleak colors
% Save data:
folder = './Subroutines/matchmaker/data';
if ~exist(folder,'dir'); mkdir(folder); end
filename = [folder '/data' datafilename '.mat'];
save(filename,'data','depth','depth_no','species','colours')
end