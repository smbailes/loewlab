%DBSCAN for volunteers 
clear all;
close all;
tic
%% Call to ChangeOverTime
Changeovertime;

%% Patient Selection
%     [location, ptID] = pathfinder; 

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end

%% Get statistics 

[r c] = size(I_mat{1}); 
for i = 1:15
    I1 = I_mat{i}(find(I_mat{i}>0));
    
    highcol = max(I1);
    high(i) = max(highcol);
    
    lowcol = min(I1);
    low(i) = min(lowcol);
    
    range(i) = high(i) - low(i);
    
    average(i) = mean2(I1);
    stdev(i) = std2(I1);
   
end     
%% Parameters
    
minPts = 10; percent = 80; 
if stdev(1) < 225
    epsilon = 6;
elseif(stdev(1) < 300 && stdev(1) >= 225)
    epsilon = 6.25;
elseif stdev(1) >= 300 
    epsilon = 6.5;
end 

if range(1) > 3000
    s = sqrt(4/3);
elseif (range(1) <= 3000 && range(1) > 2700)
    s = sqrt(10/8.5);
elseif range(1) <= 2700
    s = 1;
end     
% prompt = {'Epsilon Value Measures Cluster Closeness. Enter Epsilon Value:',...
% 'Enter MinPts:','Enter Desired %:','Enter desired scaling factor'};  
% dlg_title = 'DBSCAN Parameters';                                         % box title
% num_lines = 1;                                                          % lines per answer
% defaultans = {'5.7','12','80','sqrt(4/3)'};          % default inputs
% options.Resize = 'on';                                                  % allows for resizing of box
% answer = inputdlg(prompt, dlg_title, [1 50], defaultans, options);      % creates box
% epsilon = str2double(answer{1});                
% minPts = str2double(answer{2});                 
% percent = str2num(answer{3});
% s = str2num(answer{4});
scaling = 1/(s^2);
fprintf('Epsilon: %d \nminPts: %d \nScaling Factor: %d\n', epsilon, minPts,scaling);

%% Plot Image with Clusters using DBSCAN
%     [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID);
for n = 7:7                  % Iterate through cell matrix for each minute
    I = I_mat{n};               % Get Image
    [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID, s, percent);


    hold off;
%     ClusterInfo CELL ARRAY
    ClusterInfo{n,1} = ClustStruct;       %Cell 1 is ClusterStructure
    ClusterInfo{n,2} = I;                   %Cell2 is Image
    ClusterInfo{n,3} = ClustData;         %Cell 3 is the ClusterData output from DBSCAN
    
end    

fprintf('Finished Clustering\n');

%% Cluster Checks
for c = 7:7
   thisImage = ClusterInfo{c,1};
   pic = ClusterInfo{c,2}; %pic = I
   
   numClust = length(thisImage);
   
   %Remove bottom border
   for i = 1:numClust %Iterate through Clusters
       clustPoints = thisImage(i).ClusterIndices; %Get cluster indices
       for a = 1:length(clustPoints(:,1)) %Search through cluster indices
           if (pic((clustPoints(a,2)+1), clustPoints(a,1)) == 0) %If pixel below any cluster has intensity 0, mark cluster for removal
               thisImage(i).RemoveCluster = 1;
           end
       end       
   end    
   fprintf('Removed clusters from bottom border\n');
   
   %Remove small and large clusters
   for p = 1:numClust
      clustPoints = thisImage(p).ClusterIndices;
      for b = 1:length(clustPoints(:,1))
          if length(clustPoints(:,1)) < minPts || length(clustPoints(:,1)) > 200
              thisImage(p).RemoveCluster = 1;
          end
      end
   end
   fprintf('Removed small/large clusters\n');
   
   %Check for vessels
   for t = 1:numClust
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
    fprintf('Removed vessels\n');
    
    ClusterInfo{c,1} = thisImage; %Save Info to ClusterInfo
end  

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
    title(sprintf('%s - Post Symmetrical Cluster Analysis',ptID));
    hold on 
     
    ClusterInfo{c,3} = clustData; %Save updated Cluster Info to Array
end

%% Idenfity left or right from Changeovertime data

lowchange = mean2(lowsquarechange);
lowsquarechange = lowavesquarechange(topx);


%% Look at clusters over time

for o = 7:7
    thisImage = ClusterInfo{o,1}; %Get Current Image Info
    
    numClusters = length(thisImage);
    
    for p = 1:numClusters
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
    a = 0; b=0;c=0;
    for l = 1:numClusters
        if(abs(totalChange(l)) > abs(lowchange)) %If the total change of a cluster is too high
            thisImage(l).RemoveCluster = 1;
            a = a+1;
        end
        if(abs(avgStepChange(l)) > abs(lowsquarechange)) %If the average change is too high
            thisImage(l).RemoveCluster = 1;
            b = b+1;
        end
        for n = 1:14 %If any of the differences between 2 times is too high
            if(abs(stepChange(n,l)) > abs(lowchange))
                thisImage(l).RemoveCluster = 1;
                c = c+1;
            end
        end
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
  
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end

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
        
%         OldXWidth(t) = max(cell2mat(XDistances));
%         OldYWidth(t) = max(cell2mat(YDistances));
        CenX = ClusterInfo{7,1}(t).ClusterCentroid(1);
        CenY = ClusterInfo{7,1}(t).ClusterCentroid(2);
        counter = 1;
        
        for i = 1:length(thisImage(t).ClusterIndices)
            if ClusterInfo{7,1}(t).ClusterIndices(i,1) == CenX
                ValuesWithCenX(counter) = ClusterInfo{7,1}(t).ClusterIndices(i,2);
                counter = counter + 1;
            end
        end
        
        XWidth(t) = length(ValuesWithCenX);
        counter = 1;
        
        for i = 1:length(thisImage(t).ClusterIndices)
            if ClusterInfo{7,1}(t).ClusterIndices(i,2) == CenY
                ValuesWithCenY(counter) = ClusterInfo{7,1}(t).ClusterIndices(i,1);
                counter = counter + 1;
            end
        end
        YWidth(t) = length(ValuesWithCenY);

        [Ydimen,Xdimen] = size(I_mat{15});
        
        for j = 1:length(ClusterInfo{7,1}(t).ClusterIndices(:,1))
            NewVesClustXpos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) + (XWidth(t)+5);
            NewVesClustYpos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) + (YWidth(t)+5); 
            if NewVesClustXpos{j} > Xdimen
                NewVesClustXpos{j} = Xdimen;
            end
            if  NewVesClustYpos{j} > Ydimen
                NewVesClustYpos{j} = Ydimen;
            end
        end
        
        for j = 1:length(ClusterInfo{7,1}(t).ClusterIndices(:,1))
            NewVesClustXneg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) - (XWidth(t)+5);
            NewVesClustYneg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) - (YWidth(t)+5);  
            if NewVesClustXneg{j} < 1
                NewVesClustXneg{j} = 1;
            end
            if NewVesClustYneg{j} < 1
                NewVesClustYneg{j} = 1;
            end
        end
        
        for k = 1:length(ClusterInfo{7,1}(t).ClusterIndices(:,1))
            OriginalClusterYIndices{k} = ClusterInfo{7,1}(t).ClusterIndices(k,2);
            OriginalClusterXIndices{k} = ClusterInfo{7,1}(t).ClusterIndices(k,1);
        end


        NewVesClustIndicesPos{t} = transpose([NewVesClustXpos;NewVesClustYpos]);
        NewVesClustIndicesNeg{t} = transpose([NewVesClustXneg;NewVesClustYneg]);
        NewVesClustIndicesRight{t} = transpose([NewVesClustXpos;OriginalClusterYIndices]);
        NewVesClustIndicesLeft{t} = transpose([NewVesClustXneg;OriginalClusterYIndices]);
        
        
        clear NewVesClustXpos NewVesClustYpos NewVesClustXneg NewVesClustYneg OriginalClusterYIndices OriginalClusterXIndices
        
        AdjustedVesselsPos = I_mat{7}(cell2mat(NewVesClustIndicesPos{t}(:,2)), cell2mat(NewVesClustIndicesPos{t}(:,1)));
        AdjustedVesselsNeg = I_mat{7}(cell2mat(NewVesClustIndicesNeg{t}(:,2)), cell2mat(NewVesClustIndicesNeg{t}(:,1)));
        AdjustedVesselsRight = I_mat{7}(cell2mat(NewVesClustIndicesRight{t}(:,2)), cell2mat(NewVesClustIndicesRight{t}(:,1)));
        AdjustedVesselsLeft = I_mat{7}(cell2mat(NewVesClustIndicesLeft{t}(:,2)), cell2mat(NewVesClustIndicesLeft{t}(:,1)));
        
         AdjVesPos = AdjustedVesselsPos(find(AdjustedVesselsPos>0));
        AdjVesNeg = AdjustedVesselsNeg(find(AdjustedVesselsNeg>0));
        AdjVesRight = AdjustedVesselsRight(find(AdjustedVesselsRight>0));
        AdjVesLeft = AdjustedVesselsLeft(find(AdjustedVesselsLeft>0));
        
        
        avgAdjustedVesselNeg(t) = mean2(AdjVesPos);
        avgAdjustedVesselPos(t) = mean2(AdjVesNeg);
        avgAdjustedVesselRight(t) = mean2(AdjVesRight);
        avgAdjustedVesselLeft(t) = mean2(AdjVesLeft);
        
        VesselPDiff = avgs(7,t)*0.25;
        if((avgAdjustedVesselPos(t) + VesselPDiff < avgs(7,t)) || (avgAdjustedVesselNeg(t)+ VesselPDiff < avgs(7,t)))
            thisImage(t).RemoveCluster = 1;
        end
%         if(avgAdjustedVesselRight(t) + VesselPDiff < avgs(7,t) || avgAdjustedVesselLeft(t) + VesselPDiff < avgs(7,t))
%             thisImage(t).RemoveCluster = 1;
%         end 
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
    
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end

%% Corresponding Nipple check: Get coordinates of nipples
figure('Name','Select nipple (w/o Tumor)'), 
 for i = 1:15
    I1 = I_mat{i}(find(I_mat{i}>0));
    imshow(I_mat{i}, [min(I1) max(I1)])
    hold on
    
    [X_corrNip{i},Y_corrNip{i}] = ginput(1)
 end
close
questdlg('Switch sides','Switch sides','Ok','Sure','Ok')

figure('Name','Select nipple (w/ Tumor)'), 
for i =1:15
    I1 = I_mat{i}(find(I_mat{i}>0));
    imshow(I_mat{i}, [min(I1) max(I1)]) % gets coordinates of nipple
    hold on,
    [X_tumNip{i},Y_tumNip{i}] = ginput(1)
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
    Indices{i} = JustClust{i}(2,1)
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

AdjustedClustStruct = struct('ClusterNumber',d,'TumorBreastClustorXPoints',XClusterIndices,'TumorBreastClustorYPoints',YClusterIndices)
NewTumXPoints = cell(15,NumClust)
NewTumYPoints = cell(15,NumClust)
[rtum,ctum] = size(I_mat{7});

for i = 1:NumClust
    for j = 1:15
        NewTumXPoints{j,i} = XClusterIndices{i} + cell2mat(Xtumchange(j)); % adjusts indices by the change in the nipple relative to time 7
        NewTumYPoints{j,i} = YClusterIndices{i} + cell2mat(Ytumchange(j));
    end % New Xpoints. i = number cluster and j = points at time
end
[r,c] = size(NewTumXPoints)
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

NumberOfRemoval = [];
for z = 1:numClustersLeft
    %If the cluster changes more than the corresponding region or it is an
    %outlier
    if(ClusterDifference(z) > 0 || abs(ClusterDifference(z)) > 3000 || clusterDifferenceData(z)>0)
        counter = counter+1;
        NumberOfRemoval(counter) = z;
    end
end 


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
    
    figure('Name','Remaining Clusters'), imshow(picture,[min(pic_adj),max(pic_adj)]); %Show Image
    hold on
    PlotClusterinResult(hrEnd,clEnd); %Display remaining clusters
    title(sprintf('%s - Final Clusters',ptID));
    hold on 
    
    %%PLOT MIDLINE IMAGE

%     figure, imshow(picture, [min(pic_adj), max(pic_adj)]);
%     hold on
%     y1=get(gca,'ylim');
%     plot([xcent xcent],y1);
%         
%     title(sprintf('%s - Mirror Midline Isolation',ptID));
%     
    ClusterInfo{c,3} = clustData; %Save updated Cluster Info to Array
end
fprintf('DBSCAN_volunteers took %04f seconds to run\n',toc)
