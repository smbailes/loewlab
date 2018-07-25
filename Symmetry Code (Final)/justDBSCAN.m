clear all;


%% Patient Selection
    [location, ptID] = pathfinder; 

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=0;
    for i=1:15          
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
    notes = txt(index,7);                   % Get any notes from txt file
    celldisp(notes);   
    
%% DBSCAN Parameters
    
prompt = {'Epsilon Value Measures Cluster Closeness. Enter Epsilon Value:',...
    'Enter MinPts:','Enter Desired %:','Enter desired scaling factor'};  
dlg_title = 'DBSCAN Parameters';                                         % box title
num_lines = 1;                                                          % lines per answer
defaultans = {'6.25','10','80','sqrt(4/3)'};          % default inputs
options.Resize = 'on';                                                  % allows for resizing of box
answer = inputdlg(prompt, dlg_title, [1 50], defaultans, options);      % creates box
epsilon = str2double(answer{1});                
minPts = str2double(answer{2});                 
percent = str2num(answer{3});
s = str2num(answer{4});
scaling = 1/(s^2);
fprintf('Epsilon: %d \nminPts: %d \nScaling Factor: %d\n', epsilon, minPts,scaling);

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
      
%% ROI Identification on First Image (Remove later)

    I1 = I_mat{7};              % Display first image
    I = getMatrixOutliers(I1);  % Remove outliers
    I_adj1 = I1(find(I1>0));    % Remove zero pixels
    I_sort1 = sort(I_adj1);
    figure('Name','Nipple Identification')
    imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast

    if strcmp(sideString,'Left') == 1                   % Direct user to tumor side
        xlabel('-->')
    else
        xlabel('<--')
    end 
    hold on
    fprintf('Select Reference Nipple \n');  % User input of nipple region
    [Xorg,Yorg] = ginput(1);
    plot(Xorg,Yorg, '*');                   % Plot center of nipple
    Xnew = Xorg + dx;                       % Coordinates of center of tumor
    Ynew = Yorg - dy;    
    xbox = xbox * scale;                    % Convert X and Y dimensions of tumor from cm to pixels
    ybox = ybox * scale;                    % xbox and ybox are length of x and y sides in pixels
 
    c1 = [round(Xnew - xbox/2), round(Ynew - ybox/2)];  % Top Left Corner
    c2 = [round(Xnew - xbox/2), round(Ynew + ybox/2)];  % Bottom Left Corner
    c3 = [round(Xnew + xbox/2), round(Ynew + ybox/2)];  % Bottum Right Corner
    c4 = [round(Xnew + xbox/2), round(Ynew - ybox/2)];  % Top Right Corner      
    hold off  
    close    
%}


%% Plot Image with Clusters using DBSCAN
%     [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID);
for n = 7:7                  % Iterate through cell matrix for each minute
    I = I_mat{n};               % Get Image
    [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID, s, percent);
    plot([c1(1 ) c2(1)],[c1(2) c2(2)],'b');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'b');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'b');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'b');
    hold off;
%     ClusterInfo CELL ARRAY
    ClusterInfo{n,1} = ClustStruct;       %Cell 1 is ClusterStructure
    ClusterInfo{n,2} = I;                   %Cell2 is Image
    ClusterInfo{n,3} = ClustData;         %Cell 3 is the ClusterData output from DBSCAN
    
end    

fprintf('Finished Clustering\n');