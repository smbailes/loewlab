%% DBSCAN before threshold
%input = [xcoord, ycoord, intensity]

clear all;
close all;
clc;
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
    'Enter MinPts:','Enter Desired %:'};  
dlg_title = 'DBSCAN Parameters';                                         % box title
num_lines = 1;                                                          % lines per answer
defaultans = {'5','10','5'};          % default inputs
options.Resize = 'on';                                                  % allows for resizing of box
answer = inputdlg(prompt, dlg_title, [1 50], defaultans, options);      % creates box
epsilon = str2double(answer{1});                
minPts = str2double(answer{2});                 
percent = str2num(answer{3});            
fprintf('Epsilon: %d \nminPts: %d \n', epsilon, minPts);

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

    I1 = I_mat{1};              % Display first image
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
%% Store x,y,intensity information for all images
for n = 1:14                    % Iterate through cell matrix for each minute
    I_n = I_mat{n};               % Get Image
    I = getMatrixOutliers(I_n);   % Remove Outliers
    I_adj = I(find(I>0));
    
    [r b] = find(I);

    % Find xyz info
%     for a = 1:length(r)
%         xyz(a,1) = b(a); %X-coordinate
%         xyz(a,2) = r(a); %Y-coordinate
%         xyz(a,3) = I(r(a), b(a)); %value at that x,y coordinate
%     end  
end 
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

%{
    %% Identify vertical centerline (BY FINDING CENTER OF CROP REGION)
    % USES THE SHORTEST COLUMN OF NONZERO PIXELS as MIDLINE COLUMN
    
%     for b = 1:14
        b = 14;
        currentImage = ClusterInfo{b,2};    %Current image
        [r,c] = size(currentImage);         %Get Size of Current Image
        
        %Search Range for Midline is Middle 1/2 of the Image (Columnwise)
        boundL = floor(c / 4); 
        boundU = ceil((3*c) / 4);
        
        
        colStart = currentImage(:,boundL); %Start Search Column from CurrentImage
        colEnd = currentImage(:,boundU); %End Search Column from Current Image
        
        holdVal = numel(colStart(find(colStart > 0))); %Holdval is the number of non-zero indices in the First column
        
        for i = boundL:boundU % Check only in the middle region of the image for the smallest non-zero column
            column = currentImage(:,i);
            
            check = numel(column(find(column > 0))); %Check is number of non-zero pixels in current column
            count(i) = check;
            
            if check < holdVal %hold is the column with the fewest nonzero pixels
                holdVal = check; %If current column has fewer nonero pixels, replace holdval with current column number
            end            
            
            xcent = floor(median(find(count == holdVal))); %If numerous columns have same number of minimal pixels, take middle column
            
        end
        
        if ((xcent >= (boundU - 10)) || (xcent <=  (boundL + 10))) % if Center is found at borders of search range, reset image center to middle column
           xcent = round(c / 2); 
        end
        
        ClusterInfo{b,4} = xcent;  %Save Centerline info to Cell Matrix
%     end      
    
  
    
    %% First Cluster Check: Midline
    for c = 1:14
        thisImage = ClusterInfo{c,1}; %Current Image Struct
        
        numClusters = length(thisImage); %Number of Clusters in Image
        
        xcent = ClusterInfo{c,4}; %Obtain Centerline for current image
        
        for a = 1:numClusters
            for i = 1:length(thisImage(a).ClusterIndices) %Sort through each Index for each Cluster
                if (thisImage(a).ClusterIndices(i,1))
                    thisImage(a).RemoveCluster = 1; %Remove Cluster from Cluster Struct
                    break;
                end
            end
        end         
    ClusterInfo{c,1} = thisImage;    %Save Updated info to ClusterInfo
    
    end
%}
%% Check: Clusters on Bottom Border
% for c = 1:14
    c = 1;
   thisImage = ClusterInfo{c,1};
   pic = ClusterInfo{c,2}; %pic = I
   
   numClust = length(thisImage);
   
   for i = 1:numClust %Iterate through Clusters
       clustPoints = thisImage(i).ClusterIndices; %Get cluster indices
       for a = 1:length(clustPoints(:,1)) %Search through cluster indices
           if (pic((clustPoints(a,2)+1), clustPoints(a,1)) == 0) %If pixel below any cluster has intensity 0, mark cluster for removal
               thisImage(i).RemoveCluster = 1;
               break
           end
       end       
   end    
   
   ClusterInfo{c,1} = thisImage; %Save Info to ClusterInfo
% end  
%% Extract Cluster Data
% for c = 1:14
    thisImage = ClusterInfo{c,1}; %Get Current Image Info
    
    picture = ClusterInfo{c,2};
    pic_adj = picture(find(picture>0));
    
    clustData = ClusterInfo{c,3}; %Copy Cluster Data 
        
    numClusters = length(thisImage);
    
    for i = 1:numClusters %Sort through clusters for image
        if thisImage(i).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            thisImage(i).ClusterMeanIntensity = 0;
        end     
    end    
    
    for a = 1:numClusters
        CI(a) = thisImage(a).ClusterMeanIntensity;
    end
    
    CI_nonzero = CI(find(CI>0));
    CI_sorted = sort(CI_nonzero);
    percent1 = percent/100;
    percent_ind1 = round(percent1*numel(CI_sorted));
    percent_val = CI_sorted(end-percent_ind1);
    
    for k = 1:numClusters
        if thisImage(k).ClusterMeanIntensity < percent_val
            thisImage(k).RemoveCluster = 1;
        end
    end
   
    
    for m = 1:numClusters %Sort through clusters for image
        if thisImage(m).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            rows = find(clustData(:,3) == m);
            clustData(rows,:) = [];
        end     
    end 
    
    hrEnd = clustData(:,(1:2)); %Cluster Indices of remaining Clusters
    clEnd = clustData(:,3); %Cluster Number for remaining clusters
    
    figure, imshow(picture,[min(pic_adj),max(pic_adj)]); %Show Image
    hold on
    PlotClusterinResult(hrEnd,clEnd); %Display remaining clusters
    title(sprintf('%s - Post Symmetrical Cluster Analysis',ptID));
    plot([c1(1) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold off
    
    %%PLOT MIDLINE IMAGE

%     figure, imshow(picture, [min(pic_adj), max(pic_adj)]);
%     hold on
%     y1=get(gca,'ylim');
%     plot([xcent xcent],y1);
%     
%     plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
%     plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
%     plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
%     plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
%         
%     title(sprintf('%s - Mirror Midline Isolation',ptID));
%     
    ClusterInfo{c,3} = clustData; %Save updated Cluster Info to Array
% end



    