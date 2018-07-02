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
defaultans = {'5.2','10','80','sqrt(5/3)'};          % default inputs
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
 
    th = 0:pi/50:2*pi;
    if xbox <= ybox
        xunit = xbox * cos(th) + Xnew;
        yunit = xbox * sin(th) + Ynew;
    elseif ybox < xbox
        xunit = ybox * cos(th) + Xnew;
        yunit = ybox * sin(th) + Ynew;      
    end 
    hold off  
    close    
%}


%% Plot Image with Clusters using DBSCAN
%     [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID);
for n = 7:7                  % Iterate through cell matrix for each minute
    I = I_mat{n};               % Get Image
    [ClustStruct, ClustData] = symmetry_cluster1(I, epsilon, minPts, ptID, s, percent);
    plot(xunit, yunit);

    hold off;
%     ClusterInfo CELL ARRAY
    ClusterInfo{n,1} = ClustStruct;       %Cell 1 is ClusterStructure
    ClusterInfo{n,2} = I;                   %Cell2 is Image
    ClusterInfo{n,3} = ClustData;         %Cell 3 is the ClusterData output from DBSCAN
    
end    

fprintf('Finished Plotting Clusters\n');
%% Check: Clusters on Bottom Border
for c = 7:7
   thisImage = ClusterInfo{c,1};
   pic = ClusterInfo{c,2}; %pic = I
   
   numClust = length(thisImage);
   bb = 0;
   
   for i = 1:numClust %Iterate through Clusters
       clustPoints = thisImage(i).ClusterIndices; %Get cluster indices
       for a = 1:length(clustPoints(:,1)) %Search through cluster indices
           if (pic((clustPoints(a,2)+1), clustPoints(a,1)) == 0) %If pixel below any cluster has intensity 0, mark cluster for removal
               thisImage(i).RemoveCluster = 1;
               bb = bb+1;
           end
       end       
   end    
   
   ClusterInfo{c,1} = thisImage; %Save Info to ClusterInfo
end  
fprintf('Finished Removing Bottom Borders\n');
%% Remove small and large clusters
for d = 7:7

   thisImage = ClusterInfo{d,1};
   pic = ClusterInfo{d,2}; %pic = I
   
   numClust = length(thisImage);
   sl = 0;
   for p = 1:numClust
      clustPoints = thisImage(p).ClusterIndices;
      for b = 1:length(clustPoints(:,1))
          if length(clustPoints(:,1)) < 10 || length(clustPoints(:,1)) > 150
              thisImage(p).RemoveCluster = 1;
              sl = sl+1;
          end
      end
   end
   ClusterInfo{d,1} = thisImage;
end 
fprintf('Finished Removing Small and Large Clusters\n');

%% Check for vessels
for s = 7:7

    thisImage = ClusterInfo{s,1};
    numClusters = length(thisImage);
    vess = 0;
    for t = 1:numClusters
        indices = thisImage(t).ClusterIndices;
        xind = indices(:,2);
        yind = indices(:,1);
        
        ymax = max(yind);
        ymin = min(yind);
        ylength = ymax-ymin;
        
        xmax = max(xind);
        xmin = min(xind);
        xlength = xmax-xmin;
        
        ylengths(t,2) = ylength;
        ylengths(t,1) = t;
        if ylength > 45
            thisImage(t).RemoveCluster = 1;
            vess = vess+1;
        end
    end
    
    ClusterInfo{s,1} = thisImage;
end
fprintf('Finished Removing Vessels (by length)\n');
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
    title(sprintf('%s - Post Symmetrical Cluster Analysis',ptID));
    plot(xunit, yunit);
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

%% Select cluster to plot on histogram

% e = imfreehand(); 
% xy = wait(e); %Double click to select freehand region
% binaryImage = e.createMask(); 
% BW = uint16(binaryImage);
% figure('Name', 'Histogram with ROI');
% for n = 7:7
%     I1 = I_mat{n};
%     I2 = I_mat{n}(find(I_mat{n}>0));
% 
%     I3 = I1.*BW; %sets all pixels outside of ROI to 0 
%     I4 = I3(find(I3>0));
% 
% %     subplot(4,4,n)
%     histogram(I2,500,'FaceColor','r','EdgeColor','r');
%     hold on
%     histogram(I4,500,'FaceColor','k','EdgeColor','k');
% end


%% Idenfity left or right from Changeovertime data

if abs(totLbreastchange) < abs(totRbreastchange) &&  abs(aveLbreastchange) < abs(aveRbreastchange)
    fprintf('%s Tumor on Left\n',ptID);
    tumorSide = 'Left';
elseif abs(totLbreastchange) > abs(totRbreastchange) &&  abs(aveLbreastchange) > abs(aveRbreastchange)
    fprintf('%s Tumor on Right\n',ptID);
    tumorSide = 'Right';
else
    fprintf('Unsure\n');
end
fprintf('Tumor Truth Data: %s\n', sideString{1});

lowchange = mean2(lowsquarechange);


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

        for m = 1:14
            stepChange(m,p) = avgs(m+1,p) - avgs(m,p);
        end 
        
        avgStepChange(:,p) = mean2(stepChange(:,p));
        
    end
    
    for l = 1:numClusters
        if abs(totalChange(:,l)) > abs(lowchange)
            thisImage(l).RemoveCluster = 1;
        end
    end
    
    
%     if strcmp(tumorSide, 'Left') == 1
%         for l = 1:numClusters
%             if(abs(totalChange(:,l)) > abs(totLbreastchange) || abs(avgStepChange(:,l)) > abs(aveLbreastchange))
%                 thisImage(l).RemoveCluster = 1;
%             end 
%         end 
%     else
%         for l = 1:numClusters
%             if(abs(totalChange(:,l)) > abs(totRbreastchange) || abs(avgStepChange(:,l)) > abs(aveRbreastchange))
%                 thisImage(l).RemoveCluster = 1;
%             end 
%         end
%     end 
    
    ClusterInfo{o,1} = thisImage;
end 

%{
%         ImBinMask = zeros(size(I_mat{o}));
%         for b = 1:length(clusterIndices)
%             ImBinMask(clusterIndices(b,2),clusterIndices(b,1)) == 1;
%         end
%         BW1 = uint16(ImBinMask);
        
%         for q = 1:15
%             I1 = I_mat{q};
%             I2 = I_mat{q}(find(I_mat{q}>0));
% 
%             I3 = I1.*BW1; %sets all pixels outside of ROI to 0 
%             I4 = I3(find(I3>0));
%         end 
%}

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
    plot(xunit, yunit);
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

%}
