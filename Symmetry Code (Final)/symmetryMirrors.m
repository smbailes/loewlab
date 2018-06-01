function symmetryMirrors
%{

removes clusters that are symmetric
Created by Shannon 8/7/2017 (modification of hotRegionDBSCAN)

Last edited by Shannon: 8/7/2017 3:38p
Last edited by Aidan: 8/17/2017 5:30p


INPUTS:
- Patient ID
- Tumor Position (Side, Clock Hour, Distance, Dimensions)
- Clustering D ata (Epsilon, MinPts, Desired Intensity Percentage)

Outputs:
- Patient Image with DBSCAN Clustering 
    - Clusters Removed based on Symmetry:
        - Clusters that cross midline are removed
        - Clusters with a different cluster reflecting across midline are
            removed
        - Clusters Along bottom line of breast crop are removed

Function Reliances:
- hotRegionFinderDBSCANmod
- DBSCAN
- PlotClusterInResult
- getMatrixOutliers

%}

%% Enter Patient Data

    prompt = {'Input Patient ID Number','Input Tumor Side','Input Clock Hour:','Input Distance(cm):',...
        'Input X Tumor Dimension(cm):','Input Y Tumor Dimension(cm):',...
        'Epsilon Value Measures Cluster Closeness. Enter Epsilon Value:',...
        'Enter MinPts:','Enter Desired %:',};  % prompts
    dlg_title = 'Input Properties';                % box title
    defaultans = {'IRST009','Left','2','2.3','2.3','2.2','5','10','3.73'};           % default inputs
    options.Resize = 'on';                              % allows for resizing of box
    answer = inputdlg(prompt, dlg_title, [1 45], defaultans, options);    % creates box
    % Extract Variable Answers
    ptID = answer{1};
    sideString = answer{2};
    hr = str2double(answer{3});
    dist = str2double(answer{4});
    xbox = str2double(answer{5});
    ybox = str2double(answer{6});
    epsilon = str2double(answer{7});                   
    minPts = str2double(answer{8});
    percent = str2num(answer{9});
    
    
    
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
    rho = dist * scale;     % Distance from origin to ROI
    [dx,dy] = pol2cart(theta, rho); % Convert tumor location as angle & dist to pixel location  
      
    %% Image Input
    
    % MAC
    %   For Full images
%            location = (['/Users/AidanMurray/Data/GWBox/GRP_Loew-Doc/NadaKamona/Clinic Patients/' ptID '/']);
    %   For Cropped Images
%           location =  (['/Users/AidanMurray/Data/GWBox/GRP_Loew-Doc/NadaKamona/Clinic Patients/Manual Crop/' ptID 'Cropped/Cropped Image/']);
    
    %   For Volunteer Images
%           location = (['/Users/AidanMurray/Data/GWBox/GRP_Loew-Doc/NadaKamona/Clinic Volunteers/Manual Crop/' ptID 'Cropped/']);
          
    % WINDOWS
    %   For cropped images
%               location = (['C:\Users\Aidan\Downloads\Manual Crop\Manual Crop\' ptID 'Cropped\Cropped Image\']);
%   For full images
           location = (['D:\BreastPatients\' ptID '\']);  
    % For Volunteer Images
    
   
%  VOLUNTEERS
 %   For cropped images
%        location = (['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Volunteers\Manual Crop\' ptID 'Cropped\']);        
    %   For full images
%         location = (['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Volunteers\' ptID '\']);

    
    % Read 15 cropped images to cell matrix I_mat
%     a=0;
%     for i=1:15          
%         I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
%         a=a+120;            % Go to next image (for cropped)
%     end

% Edit IMAGE INTAKE so that two cropped mirror images compose the input
% image with midline between mirrors
im1 = imread ([location sprintf('%s_reflection.tif',ptID)]);
I_mat{1} = rgb2gray(im1);

    %% ROI Identification on First Image (Remove later) 
    
    %**ROI IGNORED FOR NOW - UNECESSARY WHEN TESTING ACCURACY WITH
    %VOLUNTEERS**
%{
    
    % Display first image
    I1 = I_mat{1};
    I = getMatrixOutliers(I1);  % Remove outliers
    I_adj1 = I1(find(I1>0));    % Non-zero pixels
    percent1 = percent/100;
    I_sort1 = sort(I_adj1);
    figure, imshow(I1,[min(I_adj1) max(I_adj1)]);
    percent_ind1 = round(percent1 * numel(I_sort1));
    percent_val1 = I_sort1(end - percent_ind1);
    [overlay_r,overlay_c] = find(I >= percent_val1); %Get Pixel locations above Percent
    if strcmp(sideString,'Left') == 1
        xlabel('-->')
    else
        xlabel('<--')
    end 
    hold on
    % Nipple Identification
    fprintf('Select Reference Nipple \n'); % User input of nipple region
    [Xorg,Yorg] = ginput(1);
    plot(Xorg,Yorg, '*');       % Plot center of nipple
    % Find Tumor Center    
    Xnew = Xorg + dx;   % Coordinates of center of tumor
    Ynew = Yorg - dy;    
    xbox = xbox * scale;    % Convert X and Y dimensions of tumor from cm to pixels
    ybox = ybox * scale;    % xbox and ybox are length of x and y sides in pixels
    c1 = [round(Xnew - xbox/2), round(Ynew - ybox/2)];  % Top Left Corner
    c2 = [round(Xnew - xbox/2), round(Ynew + ybox/2)];  % Bottom Left Corner
    c3 = [round(Xnew + xbox/2), round(Ynew + ybox/2)];  % Bottum Right Corner
    c4 = [round(Xnew + xbox/2), round(Ynew - ybox/2)];  % Top Right Corner      
    hold off
%}
    
    %% Plot Image with Clusters using DBSCAN
    
    % Won't work to run DBSCANmod: need to run straight DBSCAN
    % Copied and pasted from DBSCANmod
    
    for n = 1:1 %For Each Image - Keep at 1 for now to save on processing power
        % Get Image & Adjust Values
        I = I_mat{n};               % Get Image
%         I = getMatrixOutliers(I);   % Remove Outliers
        I_adj = I(find(I>0));       % Remove Zero Pixels
        percent = percent / 100;
        [b, edge] = histcounts(I_adj); % Get Image Histogram Data
        % Display Image Histogram w/ Percent Area
        I_sort = sort(I_adj);       % Arrange Image Hist in Order Low -> High
        percent_ind = round(percent * numel(I_sort));   % Find the index number for the User Input Percentage
        percent_val = I_sort(end - percent_ind);        % Find Intensity for the Percentage Number
        [overlay_r,overlay_c] = find(I >= percent_val); % Get Pixel locations above Percent Indicated
        % Find HotRegions
        for a = 1:length(overlay_r)
            hotregion(a,1) = overlay_c(a);
            hotregion(a,2) = overlay_r(a);
        end  
        [Clusters, isNoise] = DBSCAN(hotregion,epsilon,minPts); % Run DBSCAN on Pixels above Intensity Percentage
        ClustersNew = Clusters;
        hotregionNew = hotregion;
        ClustersNew(isNoise) = [];      % Remove all Noise Pixels from Clusters
        hotregionNew(isNoise,:) = [];   % Remove all Noise Pixels from hotRegion
        
        figure, imshow(I,[min(I_adj) max(I_adj)]);      % Display Image w Contrast
        PlotClusterinResult(hotregionNew,ClustersNew);  % Plot Clusters on Image
        title(sprintf('%s - Pre Symmetrical Cluster Analysis',ptID))
        xlabel(sprintf('Top %.2f of Pixels',percent*100))
                
        ClusterData(:,(1:2)) = hotregionNew;    % Columns 1,2 are X,Y Indices
        ClusterData(:,3) = ClustersNew;         % Column 3 is Cluster Number
        numClusters = max(ClusterData(:,3));    % Number of Clusters in Image
        
        % Allocate Cluster Structure with Necessary Fields for Data
        % Tracking
        ClusterStruct(1:numClusters) = struct('ClusterNumber',0,'ClusterIndices',[],'ClusterMeanIntensity',[],'ClusterStdIntensity',[],'StdDivMean',[],'RemoveCluster',0);
        
        % Cluster Number
        for c = 1:numClusters   
            ClusterStruct(c).ClusterNumber = c;
        end
        
        % Cluster Indices
        for a = 1:length(ClusterData) %Sort Through ClusterData Matrix
            b = ClusterData(a,3); %B is the Cluster Number of Indices in Column A
            if isempty(ClusterStruct(b).ClusterIndices) %If ClusterStruct entry has no indices, input the first indices from the cluster
                ClusterStruct(b).ClusterIndices(1,:) = ClusterData(a,(1:2));
            else %If Cluster entry already has indices, place new ones in end+1 row
                ClusterStruct(b).ClusterIndices(end+1,:) = ClusterData(a,(1:2));
            end
        end
        % Calculate Statistics 
        for i = 1:numClusters %For each cluster, calculate the Mean and Standard Deviation for Intensity 
            ClusterStruct(i).ClusterMeanIntensity = mean2(I(ClusterStruct(i).ClusterIndices(:,2),ClusterStruct(i).ClusterIndices(:,1)));
            ClusterStruct(i).ClusterStdIntensity = std2(I(ClusterStruct(i).ClusterIndices(:,2),ClusterStruct(i).ClusterIndices(:,1)));
            
            %StandardDeviation / Mean is the indicator from Li Jiang Paper 
            ClusterStruct(i).StdDivMean = ClusterStruct(i).ClusterStdIntensity / ClusterStruct(i).ClusterMeanIntensity;
             %Matrix Referral is (Row,Column) so must reverse the indices when calling the statistics
        end
%        [Struct,Image] = hotRegionFinderDBSCANmod(I_mat{b},5);
        % Cell matrix ClusterData contains Cluster Info for each Image
        
        %ClusterInfo CELL ARRAY
        ClusterInfo{n,1} = ClusterStruct;  %Cell 1 is ClusterStructure
        ClusterInfo{n,2} = I; %Cell2 is Image
        ClusterInfo{n,3} = ClusterData; %Cell 3 is the ClusterData output from DBSCAN

    end
    
    %% Identify vertical centerline (BY FINDING CENTER OF CROP REGION)
    % USES THE SHORTEST COLUMN OF NONZERO PIXELS as MIDLINE COLUMN
    
    for b = 1:1 %Iterate through the cell array for Image Cluster Data
        currentImage = ClusterInfo{b,2}; %Current image
%         [r,c] = size(currentImage); %Get Size of Current Image
        c = size(currentImage,2);
%         boundL = floor(c / 4); %Search Range for Midline is Middle 1/2 of the Image (Columnwise)
%         boundU = ceil((3*c) / 4);
%         
%         
%         colStart = currentImage(:,boundL); %Start Search Column from CurrentImage
%         colEnd = currentImage(:,boundU); %End Search Column from Current Image
%         
% %         holdVal = numel(colStart(find(colStart > 0))); %Holdval is the number of non-zero indices in the First column
%         holdVal = colStart(end);
%         
%         
%         for i = boundL:boundU % Check only in the middle region of the image for the smallest non-zero column
%             column = currentImage(:,i);
%             
% %             check = numel(column(find(column > 0))); %Check is number of non-zero pixels in current column
%             check = column(end);
%             count(i) = check;
%             
%             if check < holdVal %hold is the column with the fewest nonzero pixels
%                 holdVal = check; %If current column has fewer nonero pixels, replace holdval with current column number
%             end            
%             
%             xcent = floor(median(find(count == holdVal))); %If numerous columns have same number of minimal pixels, take middle column
%             
%         end
%         
%         if ((xcent >= (boundU - 10)) || (xcent <=  (boundL + 10))) % if Center is found at borders of search range, reset image center to middle column
%            xcent = round(c / 2); 
%         end
        xcent = round(c/2);

        ClusterInfo{b,4} = xcent;
    end      
    
  
    
    %% First Check: See if Cluster crosses center
    
    for c = 1:1
        thisImage = ClusterInfo{c,1}; %Current Image Struct
        
        numClusters = length(thisImage); %Number of Clusters in Image
        
        xcent = ClusterInfo{c,4}; %Obtain Centerline for current image
        
        for a = 1:numClusters
            for i = 1:length(thisImage(a).ClusterIndices) %Sort through each Index for each Cluster
                if (thisImage(a).ClusterIndices(i,1) == xcent || ((thisImage(a).ClusterIndices(i,1)>(xcent-5)) && (thisImage(a).ClusterIndices(i,1) <(xcent+5))))
                    %If Cluster has Indices on midline, it crosses over
%                     thisImage(a) = [];
                    thisImage(a).RemoveCluster = 1; %Remove Cluster from Cluster Struct
                    break;
                end
            end
        end
%         
    ClusterInfo{c,1} = thisImage;    %Save Updated info to ClusterInfo
    
    end
    
    
    %% Second Check: Check for Cluster Reflections (TEST)
    
    for c = 1:1        
        
        thisImage = ClusterInfo{c,1};
        
        numClusters = length(thisImage);
        
        xcent = ClusterInfo {c,4};
        
        for a = 1:(numClusters-1) %Sort through Each Cluster (save the last one, unecessarily redundant)
                       
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
%                             holdval = 1;
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
for c = 1:1
   thisImage = ClusterInfo{c,1};
   pic = ClusterInfo{c,2};
   
   numClust = length(thisImage);
   [rows,col] = size(pic);
   
   for i = 1:numClust %Iterate through Clusters
       clustPoints = thisImage(i).ClusterIndices; %Get cluster indices
       for a = 1:length(clustPoints) %Search through cluster indices
           if (clustPoints(a,2) == 1 || clustPoints(a,2) == rows || clustPoints(a,1) == 1 || clustPoints(a,1) == col)
%            if (pic((clustPoints(a,2)+1), clustPoints(a,1)) == 0) %If pixel below any cluster has intensity 0, mark cluster for removal
               thisImage(i).RemoveCluster = 1;
               break
           end
       end       
   end    
   
   
   ClusterInfo{c,1} = thisImage; %Save Info to ClusterInfo
end
    
    
    %% Extract Cluster Data


for c = 1:1
    
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
    
    figure, imshow(picture, [min(pic_adj), max(pic_adj)]);
    hold on
    y1=get(gca,'ylim')
    plot([xcent xcent],y1);
    title(sprintf('%s - Mirror Midline Isolation',ptID));
    
end


%% Identify Clusters that are symmetric 


%% End
%% Close Function
    fprintf('Press Any Key to End Function\n'); 
    pause % Wait for UI Key
%     close all
  %  clc     % Close and  Clear Function & CMD Window    
end