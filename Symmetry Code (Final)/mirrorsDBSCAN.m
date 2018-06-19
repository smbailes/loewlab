clear all;
close all;
clc;
%% Patient Selection
    [location, ptID] = pathfinder; 
    
%% Read in Images
    left0 = imread([location '\' ptID '_Left_min0.tif']);
    right0 = imread([location '\' ptID '_Right_min0.tif']);
    left15 = imread([location '\' ptID '_Left_min15.tif']);
    right15 = imread([location '\' ptID '_Right_min15.tif']);
    
    I_mat{1} = left0; I_mat{2} = right0; 
    I_mat{3} = left15; I_mat{4} = right15;
    
%% DBSCAN Parameters
    
prompt = {'Epsilon Value Measures Cluster Closeness. Enter Epsilon Value:',...
    'Enter MinPts:','Enter Desired %:','Enter desired scaling factor'};  
dlg_title = 'DBSCAN Parameters';                                         % box title
num_lines = 1;                                                          % lines per answer
defaultans = {'5','10','5','1'};          % default inputs
options.Resize = 'on';                                                  % allows for resizing of box
answer = inputdlg(prompt, dlg_title, [1 50], defaultans, options);      % creates box
epsilon = str2double(answer{1});                
minPts = str2double(answer{2});                 
percent = str2num(answer{3});
s = str2num(answer{4});
fprintf('Epsilon: %d \nminPts: %d \nScaling Factor: %d\n', epsilon, minPts,s);


%% Plot Image with Clusters using DBSCAN
n = 1;
%     [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID);
% for n = 1:14                    % Iterate through cell matrix for each minute
    I = I_mat{n};               % Get Image
    [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID);
    
    plot([c1(1) c2(1)],[c1(2) c2(2)],'k');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'k');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'k');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'k');
    hold off;
%     ClusterInfo CELL ARRAY
    ClusterInfo{n,1} = ClustStruct;       %Cell 1 is ClusterStructure
    ClusterInfo{n,2} = I;                   %Cell2 is Image
    ClusterInfo{n,3} = ClustData;         %Cell 3 is the ClusterData output from DBSCAN
    
% end    

%% Plot Image with Clusters using DBSCAN
n = 1;
%     [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID);
% for n = 1:14                    % Iterate through cell matrix for each minute
    I = I_mat{n};               % Get Image
    [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID, s);

    hold off;
%     ClusterInfo CELL ARRAY
    ClusterInfo{n,1} = ClustStruct;       %Cell 1 is ClusterStructure
    ClusterInfo{n,2} = I;                   %Cell2 is Image
    ClusterInfo{n,3} = ClustData;         %Cell 3 is the ClusterData output from DBSCAN
    
% end    