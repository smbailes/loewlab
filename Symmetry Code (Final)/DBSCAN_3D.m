%% DBSCAN before threshold
%Function Reliances:
% - pathfinder
% - patientselect
% - volunteerselect
% - DBSCAN
% - PlotClusterInResult
% - getMatrixOutliers
% - symmetry_cluster1


clear all;
close all;
clc;
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

%% Show Histogram for whole image and Truth Region
%     imshow(I_mat{7},[min(I_adj1) max(I_adj1)]);               % Display with contrast
%     hold on;
%     plot([c1(1) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
%     plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
%     plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
%     plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
%     hold on;
%     ROI = imrect();
%     xy = wait(ROI);
%     close
%     
%     figure('Name','Histograms'), subplot(4,4,1);
%     for n = 1:14
%         I2 = I_mat{n}(find(I_mat{n}>0));
%         newCrop = imcrop(I_mat{n}, xy);
% 
%         subplot(4,4,n)
%         histogram(I2,1000,'FaceColor','r','EdgeColor','r');
%         hold on
%         histogram(newCrop,1000,'FaceColor','k','EdgeColor','k');
%     end

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


%% Check: Clusters on Bottom Border
% for c = 1:14
   c = 1;
   thisImage = ClusterInfo{c,1};
   pic = ClusterInfo{c,2}; %pic = I
   
   numClust = length(thisImage);
   bb = 0;
   
   for i = 1:numClust %Iterate through Clusters
       clustPoints = thisImage(i).ClusterIndices; %Get cluster indices
       for a = 1:length(clustPoints(:,1)) %Search through cluster indices
           if (pic((clustPoints(a,2)+3), clustPoints(a,1)) == 0) %If pixel below any cluster has intensity 0, mark cluster for removal
               thisImage(i).RemoveCluster = 1;
               bb = bb+1;
               break
           end
       end       
   end    
   
   ClusterInfo{c,1} = thisImage; %Save Info to ClusterInfo
% end  

%% Remove small and large clusters
   c = 1;
   thisImage = ClusterInfo{c,1};
   pic = ClusterInfo{c,2}; %pic = I
   
   numClust = length(thisImage);
   sl = 0;
   for p = 1:numClust
      clustPoints = thisImage(p).ClusterIndices;
      for b = 1:length(clustPoints(:,1))
          if length(clustPoints(:,1)) < minPts || length(clustPoints(:,1)) > 150
              thisImage(p).RemoveCluster = 1;
              sl = sl+1;
          end
      end
   end
   
ClusterInfo{c,1} = thisImage;

%% Remove Cluster Data
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
    
     ClusterInfo{c,1} = thisImage;
%% Keep certain %
    
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
   
 ClusterInfo{c,1} = thisImage;
%% Plot left over clusters
    
    for m = 1:numClusters %Sort through clusters for image
        if thisImage(m).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            rows = find(clustData(:,3) == m);
            clustData(rows,:) = [];
        end     
    end 
    
    hrEnd = clustData(:,(1:2)); %Cluster Indices of remaining Clusters
    clEnd = clustData(:,3); %Cluster Number for remaining clusters
    
    figure('Name','Remaining Clusters'), imshow(picture,[min(pic_adj),max(pic_adj)]); %Show Image
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





    