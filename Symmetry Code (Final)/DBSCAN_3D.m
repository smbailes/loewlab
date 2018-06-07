%% DBSCAN before threshold
%input = [xcoord, ycoord, intensity]

%% Patient Selection
    [location, ptID] = pathfinder; 

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=120;
    for i=1:14          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end

%% Import Truth Data from Excel file
    [num,txt,raw] = xlsread('TruthData');   % Put TruthData into 3 cell arrays
    index = find(ismember(txt,ptID));       % Select row where patient data is
    hr = num(index-1,1);                    % Get clock hour from num file
    dist = num(index-1,2);                  % Get distance from num file
    xbox = num(index-1,3);                  % Get x-dim from num file
    ybox = num(index-1,4);                  % Get y-dim from num file
    sideString = txt(index,2);              % Get side from txt file
    notes = txt(index,8);                   % Get any notes from txt file
    celldisp(notes);   
    
%% DBSCAN Parameters
    
prompt = {'Epsilon Value Measures Cluster Closeness. Enter Epsilon Value:',...
    'Enter MinPts:'};  
dlg_title = 'DBSCAN Parameters';                                         % box title
num_lines = 1;                                                          % lines per answer
defaultans = {'5','10'};          % default inputs
options.Resize = 'on';                                                  % allows for resizing of box
answer = inputdlg(prompt, dlg_title, [1 50], defaultans, options);      % creates box
epsilon = str2double(answer{1});                
minPts = str2double(answer{2});                 
    
    % Convert Clock Hour to Angle (in rad) 
    if hr <= 9 && hr > 3            
        hr_ang = (abs((hr - 9)) * (pi/6)) + pi;
    elseif hr <= 3 
        hr_ang = (3 - hr) * (pi/6);
    elseif hr > 9
        hr_ang = pi - ((hr - 9) * (pi/6));
    end
    theta = hr_ang;     % angle in radians
    
    % Convert from Polar to Cartesian
    scale = 15;             % Pixels / CM Scale (Can be changed to UI)
    %Display pt ruler image and use ginput to get ruler size, change scale
    %accordingly
    rho = dist * scale;     % Distance from origin to ROI
    [dx,dy] = pol2cart(theta, rho); % Convert tumor location as angle & dist to pixel location  
      
    
%% Store x,y,intensity information for all images
for n = 1:14                    % Iterate through cell matrix for each minute
    I_n = I_mat{n};               % Get Image
    I = getMatrixOutliers(I_n);   % Remove Outliers
    I_adj = I(find(I>0));
    
    [r c] = find(I);

    % Find xyz info
    for a = 1:length(r)
        xyz(a,1) = c(a); %X-coordinate
        xyz(a,2) = r(a); %Y-coordinate
        xyz(a,3) = I(r(a), c(a)); %value at that x,y coordinate
    end  

    %% Plot Image with Clusters using DBSCAN

    [Clusters, isNoise] = DBSCAN(xyz,epsilon,minPts); % Run DBSCAN on Pixels above Intensity Percentage
    ClustersNew = Clusters;
    ClustersNew(isNoise) = [];      % Remove all Noise Pixels from Clusters

    figure('Name','Pre-Symmetrical Cluster Analysis')
    imshow(I,[min(I_adj) max(I_adj)]);                          % Display Image w Contrast
    hold on;
    PlotClusterinResult(ClustersNew, I); hold on;              % Plot Clusters on Image
    % plot([c1(1) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    % plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    % plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    % plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    title(sprintf('%s - Pre Symmetrical Cluster Analysis',ptID));
    xlabel(sprintf('Top %.2f of Pixels',percent2*100));
    hold off;
    hold off;
%     [ClusterStruct, ClusterData] = symmetry_cluster1(I, epsilon, minPts, percent, ptID);
%    
%     %ClusterInfo CELL ARRAY
%     ClusterInfo{n,1} = ClusterStruct;       %Cell 1 is ClusterStructure
%     ClusterInfo{n,2} = I;                   %Cell2 is Image
%     ClusterInfo{n,3} = ClusterData;         %Cell 3 is the ClusterData output from DBSCAN
end    



    