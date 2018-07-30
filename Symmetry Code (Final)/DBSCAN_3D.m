%% DBSCAN before threshold (or with high threshold) 
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

%% Call to ChangeOverTime
Changeovertime;
close all;

%% Patient Selection
%     [location, ptID] = pathfinder; 

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
    xbox = xbox*1.2;
    ybox = ybox*1.2;
    
    c1 = [round(Xnew - xbox/2), round(Ynew - ybox/2)];  % Top Left Corner
    c2 = [round(Xnew - xbox/2), round(Ynew + ybox/2)];  % Bottom Left Corner
    c3 = [round(Xnew + xbox/2), round(Ynew + ybox/2)];  % Bottum Right Corner
    c4 = [round(Xnew + xbox/2), round(Ynew - ybox/2)];  % Top Right Corner      
%     th = 0:pi/50:2*pi;
%     if xbox <= ybox
%         xunit = xbox * cos(th) + Xnew;
%         yunit = xbox * sin(th) + Ynew;
%     elseif ybox < xbox
%         xunit = ybox * cos(th) + Xnew;
%         yunit = ybox * sin(th) + Ynew;      
%     end 
    hold off  
    close    
%}



%% Plot Image with Clusters using DBSCAN
%     [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID);
for n = 7:7                  % Iterate through cell matrix for each minute
    I = I_mat{n};               % Get Image
    [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID, s, percent);
%     plot(xunit, yunit);
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

%% Try thresholding

thisImage = ClusterInfo{7,1};
numClust = length(thisImage);

for i = 1:numClust
    clusterIntensities(i) = thisImage(i).ClusterMeanIntensity; 
end 

percent = 0.8;
percent2 = 0.005;
clusterSort = sort(clusterIntensities);       % Arrange Image Hist in Order Low -> High
percent_ind_low = round(percent * numel(clusterSort));   % Find the index number for the User Input Percentage
percent_val_low = clusterSort(end - percent_ind_low)        % Find Intensity for the Percentage Number
percent_ind_high = round(percent2 * numel(clusterSort));
percent_val_high = clusterSort(end - percent_ind_high)

for j = 1:numClust
    if(clusterIntensities(j) < percent_val_low || clusterIntensities(j) > percent_val_high)
        thisImage(j).RemoveCluster = 1;
    end 
end

ClusterInfo{7,1} = thisImage;

%% Plot left over clusters
for g = 7:7 
    
    thisImage = ClusterInfo{g,1}; %Get Current Image Info
    numClusters = length(thisImage);
    clustData = ClusterInfo{g,3}; %Copy Cluster Data 

    for m = 1:numClusters %Sort through clusters for image
        if thisImage(m).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            rows = find(clustData(:,3) == m);
            clustData(rows,:) = [];
        end     
    end 
    
    hrEnd = clustData(:,(1:2)); %Cluster Indices of remaining Clusters
    clEnd = clustData(:,3); %Cluster Number for remaining clusters
    
    picture = ClusterInfo{g,2};
    pic_adj = picture(find(picture>0));
    
    figure('Name','Remaining Clusters'), imshow(picture,[min(pic_adj),max(pic_adj)]); %Show Image
    hold on
    PlotClusterinResult(hrEnd,clEnd); %Display remaining clusters
    title(sprintf('%s - Post Threshold',ptID));
%     plot(xunit, yunit);plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    
    %%PLOT MIDLINE IMAGE

%     figure, imshow(picture, [min(pic_adj), max(pic_adj)]);
%     hold on
%     y1=get(gca,'ylim');
%     plot([xcent xcent],y1);
%         
%     title(sprintf('%s - Mirror Midline Isolation',ptID));
%     
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end

%% Cluster Checks
for c = 7:7
   thisImage = ClusterInfo{c,1};
   pic = ClusterInfo{c,2}; %pic = I
   
   numClust = length(thisImage);
   
   %Remove bottom border
   for i = 1:numClust %Iterate through Clusters
       if thisImage(i).RemoveCluster == 0
           clustPoints = thisImage(i).ClusterIndices; %Get cluster indices
           for a = 1:length(clustPoints(:,1)) %Search through cluster indices
               if (pic((clustPoints(a,2)+1), clustPoints(a,1)) == 0) %If pixel below any cluster has intensity 0, mark cluster for removal
                   thisImage(i).RemoveCluster = 1;
               end
           end       
       end 
   end    
   fprintf('Removed clusters from bottom border\n');
   
   %Remove small and large clusters
   for p = 1:numClust
       if thisImage(p).RemoveCluster == 0 
          clustPoints = thisImage(p).ClusterIndices;
          for b = 1:length(clustPoints(:,1))
              if length(clustPoints(:,1)) < minPts || length(clustPoints(:,1)) > 200
                  thisImage(p).RemoveCluster = 1;
              end
          end
       end 
   end
   fprintf('Removed small/large clusters\n');
   
   %Check for vessels
    for t = 1:numClust
        if thisImage(t).RemoveCluster == 0
            indices = thisImage(t).ClusterIndices;
            xind = indices(:,1);
            yind = indices(:,2);

            ymax = max(yind);
            ymin = min(yind);
            ylength = ymax-ymin;

            xmax = max(xind);
            xmin = min(xind);
            xlength = xmax-xmin;

            DiagnolLength = sqrt(xlength^2 + ylength^2);

            if(ylength > 30 || xlength > 30 || DiagnolLength > 20)
                thisImage(t).RemoveCluster = 1;
            end

            if(ylength >= xlength*5 || xlength>= ylength*5)
                thisImage(t).RemoveCluster = 1;
            end 
        end 
    end
    fprintf('Removed vessels\n');
    
    ClusterInfo{c,1} = thisImage; %Save Info to ClusterInfo
end  

%% Remove Cluster Data
%{
for e = 7:7
    thisImage = ClusterInfo{e,1}; %Get Current Image Info
    clustData = ClusterInfo{e,3}; %Copy Cluster Data     
    numClusters = length(thisImage);
    
    for i = 1:numClusters %Sort through clusters for image
        if thisImage(i).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            thisImage(i).ClusterMeanIntensity = 0;
        end     
    end    
     ClusterInfo{e,1} = thisImage;
end 
fprintf('Finished Removing Cluster Mean Intensity Data\n');

%% Keep certain % of clusters
% for f = 1:14:15
%     
%     thisImage = ClusterInfo{f,1}; %Get Current Image Info
%     numClusters = length(thisImage);
%     for a = 1:numClusters
%         CI(a) = thisImage(a).ClusterMeanIntensity;
%     end
%     
%     CI_nonzero = CI(find(CI>0));
%     CI_sorted = sort(CI_nonzero);
%     percent1 = percent/100;
%     percent_ind1 = round(percent1*numel(CI_sorted));
%     percent_val = CI_sorted(end-percent_ind1);
%     
%     for k = 1:numClusters
%         if thisImage(k).ClusterMeanIntensity < percent_val
%             thisImage(k).RemoveCluster = 1;
%         end
%     end
%  ClusterInfo{f,1} = thisImage;
% end 
% fprintf('Finished Removing Clusters Below Threshold\n');
%} 

%% Plot left over clusters
for g = 7:7 
    
    thisImage = ClusterInfo{g,1}; %Get Current Image Info
    numClusters = length(thisImage);
    clustData = ClusterInfo{g,3}; %Copy Cluster Data 

    for m = 1:numClusters %Sort through clusters for image
        if thisImage(m).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            rows = find(clustData(:,3) == m);
            clustData(rows,:) = [];
        end     
    end 
    
    hrEnd = clustData(:,(1:2)); %Cluster Indices of remaining Clusters
    clEnd = clustData(:,3); %Cluster Number for remaining clusters
    
    picture = ClusterInfo{g,2};
    pic_adj = picture(find(picture>0));
    
    figure('Name','Remaining Clusters'), imshow(picture,[min(pic_adj),max(pic_adj)]); %Show Image
    hold on
    PlotClusterinResult(hrEnd,clEnd); %Display remaining clusters
    title(sprintf('%s - Post Spatial Analysis',ptID));
%     plot(xunit, yunit);plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    
    %%PLOT MIDLINE IMAGE

%     figure, imshow(picture, [min(pic_adj), max(pic_adj)]);
%     hold on
%     y1=get(gca,'ylim');
%     plot([xcent xcent],y1);
%         
%     title(sprintf('%s - Mirror Midline Isolation',ptID));
%     
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end


%% Idenfity left or right from Changeovertime data

if(abs(totLbreastchange) < abs(totRbreastchange) &&  abs(aveLbreastchange) < abs(aveRbreastchange))
    fprintf('%s Tumor on Left\n',ptID);
    tumorSide = 'Left';
elseif(abs(totLbreastchange) > abs(totRbreastchange) &&  abs(aveLbreastchange) > abs(aveRbreastchange))
    fprintf('%s Tumor on Right\n',ptID);
    tumorSide = 'Right';
else
    fprintf('Unsure\n');
end
fprintf('Tumor Truth Data: %s\n', sideString{1});

lowchange = mean2(lowsquarechange);
lowsquarechange = lowavesquarechange(topx);

%% Look at clusters over time

for o = 7:7
    thisImage = ClusterInfo{o,1}; %Get Current Image Info
    
    numClusters = length(thisImage);
    
    for p = 1:numClusters
        if thisImage(p).RemoveCluster == 0
            clusterIndices = thisImage(p).ClusterIndices;
            for q = 1:15
                I1 = I_mat{q};
                %averages of same cluster over time
                avgs(q,p) = mean2(I1(clusterIndices(:,2),clusterIndices(:,1)));
            end 
           totalChange(:,p) = avgs(1,p) - avgs(15,p);

           %Amount of change each minute 
           for m = 1:14
                stepChange(m,p) = avgs(m+1,p) - avgs(m,p);
           end 

           %Average amount of change over the 15 minutes
           avgStepChange(:,p) = mean2(stepChange(:,p));
        end   
        
    end
    for l = 1:numClusters
        if thisImage(l).RemoveCluster == 0
            if(totalChange(l) > abs(lowchange)) %If the total change of a cluster is too high
                thisImage(l).RemoveCluster = 1;
            end
            if(abs(avgStepChange(l)) > abs(lowsquarechange)) %If the average change is too high
                thisImage(l).RemoveCluster = 1;
            end
        end
%         for n = 1:14 %If any of the differences between 2 times is too high
%             if(abs(stepChange(n,l)) > abs(lowchange))
%                 thisImage(l).RemoveCluster = 1;
%             end
%         end
    end
    
    ClusterInfo{o,1} = thisImage;
end 
fprintf('Removed clusters based on change over time\n');

%% Plot left over clusters (again)  
for g = 7:7 
    
    thisImage = ClusterInfo{g,1}; %Get Current Image Info
    numClusters = length(thisImage);
    clustData = ClusterInfo{g,3}; %Copy Cluster Data 

    for m = 1:numClusters %Sort through clusters for image
        if thisImage(m).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            rows = find(clustData(:,3) == m);
            clustData(rows,:) = [];
        end     
    end 
    
    hrEnd = clustData(:,(1:2)); %Cluster Indices of remaining Clusters
    clEnd = clustData(:,3); %Cluster Number for remaining clusters
    
    picture = ClusterInfo{g,2};
    pic_adj = picture(find(picture>0));
    
    figure('Name','Remaining Clusters 2.0'), imshow(picture,[min(pic_adj),max(pic_adj)]); %Show Image
    hold on
    PlotClusterinResult(hrEnd,clEnd); %Display remaining clusters
    title(sprintf('%s - Post Cluster Intensity Analysis',ptID));
%     plot(xunit, yunit);
    plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end
%}

%% Vessel Check
for o = 7:7 
    thisImage = ClusterInfo{o,1}; %Get Current Image Info
    
    numClusters = length(thisImage);
    
  counter = 0;
  for t = 1:numClust
      
      if thisImage(t).RemoveCluster == 0
        indices = thisImage(t).ClusterIndices;
        xind = indices(:,1);
        yind = indices(:,2);
        
        ymax = max(yind);
        ymin = min(yind);
        ylength = ymax-ymin;
        
        xmax = max(xind);
        xmin = min(xind);
        xlength = xmax-xmin;
        
        DiagnolLength = sqrt(xlength^2 + ylength^2);
        ClusterSlope = ylength/xlength;
        ClusterPerpSlope = -1/ClusterSlope;
        YIntercept = (-ClusterPerpSlope*thisImage(t).ClusterCentroid(1)) + thisImage(t).ClusterCentroid(2);
        counter = 1;
        for i = 1:length(thisImage(t).ClusterIndices)
            v = floor(thisImage(t).ClusterIndices(i,2) - (ClusterPerpSlope*thisImage(t).ClusterIndices(i,1) + YIntercept));
            ClusterLineIndices{i,t} = v;
            if v == 0
              XDistances{counter}= abs((ClusterInfo{7,1}(t).ClusterIndices(i,1) - ClusterInfo{7,1}(t).ClusterCentroid(1))*2);
              YDistances{counter} = abs((ClusterInfo{7,1}(t).ClusterIndices(i,2) - ClusterInfo{7,1}(t).ClusterCentroid(2))*2);
              counter = counter + 1;  
            end
        end
     XWidth(t) = max(cell2mat(XDistances));
     YWidth(t) = max(cell2mat(YDistances));
     [Xdimen,Ydimen] = size(I_mat{15});
      for j = 1:length(ClusterInfo{7,1}(t).ClusterIndices(:,1))
         NewVesClustXpos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) + (XWidth(t)+5);
         NewVesClustYpos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) + (YWidth(t)+5);    
         if NewVesClustXpos{j} > Xdimen
             NewVesClustXpos{j} = Xdimen;
         elseif  NewVesClustYpos{j} > Ydimen
             NewVesClustYpos{j} = Ydimen;
         end
      end
      for j = 1:length(ClusterInfo{7,1}(t).ClusterIndices(:,1))
         NewVesClustXneg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) - (XWidth(t)+5);
         NewVesClustYneg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) - (YWidth(t)+5); 
         if NewVesClustXneg{j} < 1
             NewVesClustXneg{j} = 1;
         elseif NewVesClustYneg{j} < 1
             NewVesClustYneg{j} = 1;
         end
      end
      
         NewVesClustIndicesPos{t} = transpose([NewVesClustXpos;NewVesClustYpos]);
         NewVesClustIndicesNeg{t} = transpose([NewVesClustXneg;NewVesClustYneg]);
         clear NewVesClustXpos NewVesClustYpos NewVesClustXneg NewVesClustYneg
      AdjustedVesselsPos = I_mat{7}(cell2mat(NewVesClustIndicesPos{t}(:,2)), cell2mat(NewVesClustIndicesPos{t}(:,1)));
      AdjustedVesselsNeg = I_mat{7}(cell2mat(NewVesClustIndicesNeg{t}(:,2)), cell2mat(NewVesClustIndicesNeg{t}(:,1)));
      avgAdjustedVesselNeg(t) = mean2(AdjustedVesselsNeg);
      avgAdjustedVesselPos(t) = mean2(AdjustedVesselsPos);
      VesselPDiff = avgs(7,t)*0.25
%       if((avgAdjustedVesselPos(t) + VesselPDiff < avgs(7,t)) || (avgAdjustedVesselNeg(t)+ VesselPDiff < avgs(7,t)))
%           thisImage(t).RemoveCluster = 1;
%       end
    end 
  end
  ClusterInfo{o,1} = thisImage;
end 

%% Plot leftover Clusters
for g = 7:7 
    
    thisImage = ClusterInfo{g,1}; %Get Current Image Info
    numClusters = length(thisImage);
    clustData = ClusterInfo{g,3}; %Copy Cluster Data 

    for m = 1:numClusters %Sort through clusters for image
        if thisImage(m).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            rows = find(clustData(:,3) == m);
            clustData(rows,:) = [];
        end     
    end 
    
    hrEnd = clustData(:,(1:2)); %Cluster Indices of remaining Clusters
    clEnd = clustData(:,3); %Cluster Number for remaining clusters
    
    picture = ClusterInfo{g,2};
    pic_adj = picture(find(picture>0));
    
    figure('Name','Remaining Clusters 3.0'), imshow(picture,[min(pic_adj),max(pic_adj)]); %Show Image
    hold on
    PlotClusterinResult(hrEnd,clEnd); %Display remaining clusters
    title(sprintf('%s - Post Vessel Check',ptID));
%     plot(xunit, yunit);
    plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end



 %% Corresponding Nipple check: Get coordinates of nipples

figure('Name','Select nipple (w/o Tumor)'), 
 for i = 1:15
    I1 = I_mat{i}(find(I_mat{i}>0));
    imshow(I_mat{i}, [min(I1) max(I1)])
    hold on
    if strcmp(sideString,'Left') == 1                   % Direct user to tumor side
        xlabel('<--')
    else
        xlabel('-->')
    end % gets coordinates of nipple
    [X_corrNip{i},Y_corrNip{i}] = ginput(1);
 end
close
questdlg('Switch sides','Switch sides','Ok','Sure','Ok')

figure('Name','Select nipple (w/ Tumor)'), 
for i =1:15
    I1 = I_mat{i}(find(I_mat{i}>0));
    imshow(I_mat{i}, [min(I1) max(I1)]) % gets coordinates of nipple
    hold on,
    if strcmp(sideString,'Left') == 1                   % Direct user to tumor side
        xlabel('-->')
    else
        xlabel('<--')
    end 
    [X_tumNip{i},Y_tumNip{i}] = ginput(1);
end
close
%% Track clusters over time
% numClust = length(ClustStruct);
clear ClustInfoCell RemoveCluster NumClust JustClust Indices XClusterIndices YClusterIndices...
    Xtumchange Ytumchange Xcorrchange Ycorrchange AdjustedClustStruct Values Xnip2tum Ynip2tum...
    XCorrIndices YCorrIndices
ClustInfoCell = struct2cell(ClusterInfo{7,1}); %converts data into cell array
RemoveCluster = cell2mat(ClustInfoCell(7,:,:)); %separates the data indicating to remove clusters
counter = 0;
for i = 1:length(ClustInfoCell)
    if RemoveCluster(i) == 0
        counter = counter+1;
        JustClust{counter} = ClustInfoCell(:,:,i);
    end
end
   NumClust = numel(JustClust); %finds number of clusters
   clear d
for i =1:length(JustClust)
    d{i} = i;
end

for i = 1:length(JustClust)
    Indices{i} = JustClust{i}(2,1);
end

for i = 1:length(JustClust)
    XClusterIndices{i} = Indices{1,i}{1,1}(:,1);
    YClusterIndices{i} = Indices{1,i}{1,1}(:,2);
end

for i = 1:15
    Xtumchange{i} = X_tumNip{i} - X_tumNip{7};
    Ytumchange{i} = Y_tumNip{i} - Y_tumNip{7};
    Xcorrchange{i} = X_corrNip{i} - X_corrNip{7};
    Ycorrchange{i} = Y_corrNip{i} - Y_corrNip{7};
end

AdjustedClustStruct = struct('ClusterNumber',d,'TumorBreastClustorXPoints',XClusterIndices,'TumorBreastClustorYPoints',YClusterIndices);
NewTumXPoints = cell(15,NumClust);
NewTumYPoints = cell(15,NumClust);
[rtum,ctum] = size(I_mat{7});

for i = 1:NumClust
    for j = 1:15
        NewTumXPoints{j,i} = XClusterIndices{i} + cell2mat(Xtumchange(j)); % adjusts indices by the change in the nipple relative to time 7
        NewTumYPoints{j,i} = YClusterIndices{i} + cell2mat(Ytumchange(j));
    end % New Xpoints. i = number cluster and j = points at time
end
[r,c] = size(NewTumXPoints);
for i = 1:c %cluster  
Xnip2tum{i} = X_tumNip{7} - XClusterIndices{i};
Ynip2tum{i} = Y_tumNip{7} - YClusterIndices{i};
end
for i = 1:NumClust
    for j = 1:15
       XCorrIndices{j,i} = X_corrNip{j} + Xnip2tum{i};
       YCorrIndices{j,i} = Y_corrNip{j} + Ynip2tum{i};
    end
end
for j = 1:c % cluster
    for i = 1:r % time
    XTumIndice = cell2mat(NewTumXPoints(i,j)); %Creates the points for the cluster at this time
    YTumIndice = cell2mat(NewTumYPoints(i,j));
    XCorrIndice = cell2mat(XCorrIndices(i,j));
    YCorrIndice = cell2mat(YCorrIndices(i,j));
    L = length(XTumIndice);
      for k = 1:L % pixel
          if XTumIndice(k) <1 % adjusts indices if exceed image
              XTumIndice(k) = 1;              
          elseif XTumIndice(k) > ctum
              XTumIndice(k) = ctum;  
          end
          if XCorrIndice(k) <1
              XCorrIndice(k) = 1;
          elseif XCorrIndice(k) > ctum
              XCorrIndice(k) = ctum;
          end
          if YTumIndice(k) <1
              YTumIndice(k) = 1;
          elseif YTumIndice(k) > rtum
              YTumIndice(k) = rtum;
          end
          if YCorrIndice(k) <1
              YCorrIndice(k) = 1;
          elseif YCorrIndice(k) > rtum
              YCorrIndice(k) = rtum;
          end
          
          TumValues{k} = I_mat{i}(floor(YTumIndice(k)),floor(XTumIndice(k))); %Records the value at each time
          CorrValues{k} = I_mat{i}(floor(YCorrIndice(k)),floor(XCorrIndice(k)));
      end
      
      TimeClusterData{i,j} = mean(cell2mat(TumValues)); % Records the average cluster value at each time
      CorrData{i,j} = mean(cell2mat(CorrValues)); %Records average cluster value for region opposite cluster at each time
      

      clear TumValues CorrValues
    end
end
      
%% Check: Compare cluster info to region on opposite breast

numClustersLeft = length(TimeClusterData(2,:));
for y = 1:numClustersLeft
    clusterDifferenceData(y) = TimeClusterData{15,y} - TimeClusterData{1,y};
    corrRegionDifference(y) = CorrData{15,y} - CorrData{1,y};
    ClusterDifference(y) = abs(clusterDifferenceData(y))-abs(corrRegionDifference(y));
end 
counter = 0;
ClusterDifference1 = ClusterDifference(find(ClusterDifference<3000));
cutoff = std2(ClusterDifference1);

for x = 1:numClustersLeft
    for v = 1:14
        stepChangeCluster(x,v) = TimeClusterData{v+1,x} - TimeClusterData{v,x};
        stepChangeCorr(x,v) = CorrData{v+1,x} - CorrData{v,x};
    end
    aveStepChangeCluster(x) = mean2(stepChangeCluster(x,:));
    aveStepChangeCorr(x) = mean2(stepChangeCorr(x,:));
    stepChangeDifference(x) = abs(aveStepChangeCluster(x)) - abs(aveStepChangeCorr(x));
end

NumberOfRemoval = [];
for z = 1:numClustersLeft
    %If the cluster changes more than the corresponding region or it is an
    %outlier
    if((ClusterDifference(z) > 75 && clusterDifferenceData(z) < 8000) || (clusterDifferenceData(z)>0 && clusterDifferenceData(z) < 8000))
        counter = counter+1;
        NumberOfRemoval(counter) = z;
    end
end 

% counter1 = 0;
% for z = 1:numClustersLeft
%     if abs(aveStepChangeCluster(z)) > abs(aveStepChangeCorr(z))
%         counter1 = counter1+1;
%         NumberOfRemoval2(counter1) = z;
%     end 
% end

%% Remove Clusters based on corresponding region change
thisImage = ClusterInfo{7,1};
ClusterCounter = 1;
RemoveClusterCounter = 1;
for i = 1:length(thisImage);
   if thisImage(i).RemoveCluster == 0
       if RemoveClusterCounter < length(NumberOfRemoval)
          if ClusterCounter == NumberOfRemoval(RemoveClusterCounter)
            RemoveClusterCounter = RemoveClusterCounter + 1;
            thisImage(i).RemoveCluster = 1; %Mark Cluster for Removal
          end
      end
       ClusterCounter = ClusterCounter+1;
   end
end

% ClusterCounter = 1;
% RemoveClusterCounter = 1;
% for i = 1:length(thisImage);
%    if thisImage(i).RemoveCluster == 0
%        if RemoveClusterCounter < length(NumberOfRemoval2)
%           if ClusterCounter == NumberOfRemoval2(RemoveClusterCounter)
%             RemoveClusterCounter = RemoveClusterCounter + 1;
%             thisImage(i).RemoveCluster = 1; %Mark Cluster for Removal
%           end
%       end
%        ClusterCounter = ClusterCounter+1;
%    end
% end



%%     
for g = 7:7 
    
%     thisImage = ClusterInfo{g,1}; %Get Current Image Info
%     numClusters = length(thisImage);
%     clustData = ClusterInfo{g,3}; %Copy Cluster Data 

    for m = 1:numClusters %Sort through clusters for image
        if thisImage(m).RemoveCluster == 1 %Remove Indices from Marked Clusters from Cluster Data copy
            rows = find(clustData(:,3) == m);
            clustData(rows,:) = [];
        end     
    end 
    
    hrEnd = clustData(:,(1:2)); %Cluster Indices of remaining Clusters
    clEnd = clustData(:,3); %Cluster Number for remaining clusters
    
    picture = ClusterInfo{g,2};
    pic_adj = picture(find(picture>0));
    
    figure('Name','Remaining Clusters 4.0'), imshow(picture,[min(pic_adj),max(pic_adj)]); %Show Image
    hold on
    PlotClusterinResult(hrEnd,clEnd); %Display remaining clusters
    title(sprintf('%s - Final Clusters',ptID));
%     plot(xunit, yunit);
    plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    
    %%PLOT MIDLINE IMAGE

%     figure, imshow(picture, [min(pic_adj), max(pic_adj)]);
%     hold on
%     y1=get(gca,'ylim');
%     plot([xcent xcent],y1);
%         
%     title(sprintf('%s - Mirror Midline Isolation',ptID));
%     
    ClusterInfo{7,1} = thisImage;
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end

