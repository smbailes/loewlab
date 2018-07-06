%DBSCAN for volunteers 


clear all;
close all;

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
    
    
prompt = {'Epsilon Value Measures Cluster Closeness. Enter Epsilon Value:',...
'Enter MinPts:','Enter Desired %:','Enter desired scaling factor'};  
dlg_title = 'DBSCAN Parameters';                                         % box title
num_lines = 1;                                                          % lines per answer
defaultans = {'5.3','12','80','sqrt(5/3)'};          % default inputs
options.Resize = 'on';                                                  % allows for resizing of box
answer = inputdlg(prompt, dlg_title, [1 50], defaultans, options);      % creates box
epsilon = str2double(answer{1});                
minPts = str2double(answer{2});                 
percent = str2num(answer{3});
s = str2num(answer{4});
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

        if(ylength > 20 || xlength > 20)
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

