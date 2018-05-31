function ClusterInfo = symmetryTEST
%{

removes clusters that are symmetric
Created by Shannon 8/7/2017 (modification of hotRegionDBSCAN)

Last edited by Aidan: 5/31/2018


INPUTS:
- Patient ID
- Tumor Position (Side, Clock Hour, Distance, Dimensions)
- Clustering Data (Epsilon, MinPts, Desired Intensity Percentage)

Outputs:
- Patient Image with DBSCAN Clustering 
    - Clusters Removed based on Symmetry:
        - Clusters that cross midline are removed
        - Clusters with a different cluster reflecting across midline are
            removed
        - Clusters Along bottom line of breast crop are removed

Function Reliances:
- userselect
- patientselect
- volunteerselect
- DBSCAN
- PlotClusterInResult
- getMatrixOutliers
- symmetry_cluster
- symmetry_centerline
- symmetry_alignClusters

%}

%% Enter Patient Data

    close all
%% Patient Selection
answer = questdlg('ID Patient Type:','Patient Type','Patient','Volunteer','Patient');

if (strcmp(answer, 'Patient'))
    ptID = patientselect;    % Dialog Box for patient selection
    %% User Selection (CHANGE THIS TO NEW USERS / DIRECTORY REF SYSTEM)
    user = userselect;          % Dialog Box for user selection
    a = strcmp(user,'Sydney');   % Compare user input to 'Sydney"
    b = strcmp(user,'Samhita');
    c = strcmp(user,'Jacob');
    if (a == 1 && b == 0)                 % If Sydney was selected
        location = (['C:\Users\smbailes\Documents\GitHub\loewlab\Symmetry Code (Final)\Images\' ptID '\Cropped Image\']);
    elseif (a == 0 && b == 1)  % If Samhita was selected
         location = (['C:\Users\samhitamurthy\GitHub\loewlab\Symmetry Code (Final)\Images\' ptID '\Cropped Image\']);
    elseif (a == 0 && b == 0 && c == 1)
        location = (['C:\Users\Jacob\Documents\GitHub\loewlab\Symmetry Code (Final)\Images\' ptID '\Cropped Image\']);
    end
else
    ptID = volunteerselect;
    user = userselect;
    
    a = strcmp(user,'Aidan');   % Compare user input to 'Aidan"
    b = strcmp(user,'Lovelace');
    if (a == 1 && b == 0)                 % If Aidan was selected
        location = (['/Users/AidanMurray/Data/GWBox/GRP_Loew-Doc/NadaKamona/Clinic Volunteers/' ptID '/']);
    elseif (a == 0 && b == 0)
        location = (['\Users\shann\Box\GRP_Loew-Doc\NadaKamona\Clinic Volunteers\Manual Crop\' ptID 'Cropped\Cropped Image\']);
    elseif (a == 0 && b == 1)
        location = (['D:\BreastPatients\Manual Crop\' ptID 'Cropped\']); %GET ACTUAL LOCATION ON LAB PC
    end
    
    
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


%% Set DBSCAN Parameters
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
      
    
    
    %% Image Input
    
    % Read 15 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end
    
%% ROI Identification on First Image (Remove later)
    I1 = I_mat{1};              % Display first image
    I = getMatrixOutliers(I1);  % Remove outliers
    I_adj1 = I1(find(I1>0));    % Remove zero pixels
    percent1 = percent/100;
    I_sort1 = sort(I_adj1);
    figure('Name','Nipple Identification')
    imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast
    percent_ind1 = round(percent1 * numel(I_sort1));
    percent_val1 = I_sort1(end - percent_ind1);
    [overlay_r,overlay_c] = find(I >= percent_val1);    % Get Pixel locations above Percent
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
    
    %% Plot Image with Clusters using DBSCAN
for n = 1:15                    % Iterate through cell matrix for each minute
    I = I_mat{n};               % Get Image
    [ClusterStruct, ClusterData] = symmetry_cluster(I, epsilon, minPts, percent, ptID);
   
    %ClusterInfo CELL ARRAY
    ClusterInfo{n,1} = ClusterStruct;       %Cell 1 is ClusterStructure
    ClusterInfo{n,2} = I;                   %Cell2 is Image
    ClusterInfo{n,3} = ClusterData;         %Cell 3 is the ClusterData output from DBSCAN
end    
    %% Identify vertical centerline (BY FINDING CENTER OF CROP REGION)
    % USES THE SHORTEST COLUMN OF NONZERO PIXELS as MIDLINE COLUMN
    
    for b = 1:15
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
    end      
    
  
    
    %% First Cluster Check: Midline
    for c = 1:15
        thisImage = ClusterInfo{c,1}; %Current Image Struct
        
        numClusters = length(thisImage); %Number of Clusters in Image
        
        xcent = ClusterInfo{c,4}; %Obtain Centerline for current image
        
        for a = 1:numClusters
            for i = 1:length(thisImage(a).ClusterIndices) %Sort through each Index for each Cluster
                if (thisImage(a).ClusterIndices(i,1) == xcent || ((thisImage(a).ClusterIndices(i,1)>(xcent-5)) && (thisImage(a).ClusterIndices(i,1) <(xcent+5))))
                    thisImage(a).RemoveCluster = 1; %Remove Cluster from Cluster Struct
                    break;
                end
            end
        end         
    ClusterInfo{c,1} = thisImage;    %Save Updated info to ClusterInfo
    
    end
    
    
    %% Second Check: Check for Cluster Reflections
    for c = 1:15    
        
        thisImage = ClusterInfo{c,1};
        
        numClusters = length(thisImage);
        
        xcent = ClusterInfo {c,4};
        
        for a = 1:(numClusters-1) %Sort through Each Cluster (minus the last one, unecessarily redundant)
                       
            for i = 1:length(thisImage(a).ClusterIndices) %Sort through entire cluster 
                
                x_pt = thisImage(a).ClusterIndices(i,1); % x_pt is current point X index
                y_pt = thisImage(a).ClusterIndices(i,2);
                refl_dist = (xcent - x_pt); %reflection distance (distance from point to midline)
                
                x_refl = x_pt + (2*refl_dist); %X index of point when reflected across midline
                nextClust = 1;
                
                
                while (a + nextClust) <= numClusters %When cluster ID + next cluster indicater is less than or equal to total number of clusters
                    for n = 1:length(thisImage(a+nextClust).ClusterIndices) %sort through each index in next cluster
                        
                        
                        %if next cluster has index that matches the reflected
                        %index of prior cluster, both Clusters are marked for
                        %removal
                        if (thisImage(a+nextClust).ClusterIndices(n,1) == x_refl) && (thisImage(a+nextClust).ClusterIndices(n,2) == y_pt)
                            thisImage(a).RemoveCluster = 1; 
                            thisImage(a+nextClust).RemoveCluster = 1;
                            break
                        end
                    end
                    nextClust = nextClust + 1; %Iterate through reference cluster
                end
            end
        end
        
    ClusterInfo{c,1} = thisImage;    %save thisImage to ClusterInfo
        
    end
    
    
    
    %% Third Check: Clusters on Bottom Border
for c = 1:15
   thisImage = ClusterInfo{c,1};
   pic = ClusterInfo{c,2};
   
   numClust = length(thisImage);
   
   for i = 1:numClust %Iterate through Clusters
       clustPoints = thisImage(i).ClusterIndices; %Get cluster indices
       for a = 1:length(clustPoints) %Search through cluster indices
           if (pic((clustPoints(a,2)+1), clustPoints(a,1)) == 0) %If pixel below any cluster has intensity 0, mark cluster for removal
               thisImage(i).RemoveCluster = 1;
               break
           end
       end       
   end    
   
   
   ClusterInfo{c,1} = thisImage; %Save Info to ClusterInfo
end
    
    
    %% Extract Cluster Data
for c = 1:15
    
    thisImage = ClusterInfo{c,1}; %Get Current Image Info
    
    picture = ClusterInfo{c,2};
    pic_adj = picture(find(picture>0));
    
    clustData = ClusterInfo{c,3}; %Copy Cluster Data 
        
    numClusters = length(thisImage);
    
    for i = 1:numClusters %Sort through clusters for image
        if thisImage(i).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            rows = find(clustData(:,3) == i);
            clustData(rows,:) = [];
        end     
    end    
    
    hrEnd = clustData(:,(1:2)); %Cluster Indices of remaining Clusters
    clEnd = clustData(:,3); %Cluster Number for remaining clusters
    
    figure, imshow(picture,[min(pic_adj),max(pic_adj)]); %Show Image
    hold on
    PlotClusterinResult(hrEnd,clEnd); %Display remaining clusters
    title(sprintf('%s - Post Symmetrical Cluster Analysis',ptID));
    plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
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
end

%% FIND CONSISTENT CLUSTERS & SHOW HEAT INFO
[newClustInfo, maxClusters] = symmetry_alignClusters(ClusterInfo);

ClusterHeats = [];

for n = 1:15
    for o = 1:length(newClustInfo{n,1})
        ClusterHeats(newClustInfo{n,1}(o).NormalizedCluster,n) = newClustInfo{n,1}(o).ClusterMeanIntensity;
    end 
end

ClusterHeats(find(ClusterHeats == 0)) = NaN;

figure
lgnd = {};
for a = 1:maxClusters
    lgnd{end+1} = ['Cluster #' num2str(a)];
    plot(ClusterHeats(a,:))
    hold on
end
legend(lgnd,'Location','NorthEastOutside')
xlabel('Minutes')
ylabel('Intensity / Temperature')
title('Cluster Intensity vs. Time')
hold off


%% End
%% Close Function
    fprintf('Press Any Key to End Function\n'); 
    pause % Wait for UI Key
%     close all
  %  clc     % Close and  Clear Function & CMD Window    
end