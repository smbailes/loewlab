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

%Start timer
tic         

%% Patient Selection
    [location, ptID] = pathfinder; 

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end

%% Get statistics for setting DBSCAN parameters

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
minPts = 10; percent = 95; 

% epsilon = 6.25;
% s = sqrt(10/8); 
if stdev(1) < 215
    epsilon = 6;
elseif(stdev(1) < 300 && stdev(1) >= 215)
    epsilon = 6.5;
elseif stdev(1) >= 300 
    epsilon = 6.75;
end 

if range(1) > 3000
    s = sqrt(4/3);
elseif (range(1) <= 3000 && range(1) > 2700)
    s = sqrt(10/8);
elseif range(1) <= 2700
    s = 1;
end     
%{
% prompt = {'Epsilon Value Measures Cluster Closeness. Enter Epsilon Value:',...
%     'Enter MinPts:','Enter Desired %:','Enter desired scaling factor'};  
% dlg_title = 'DBSCAN Parameters';                                         % box title
% num_lines = 1;                                                          % lines per answer
% defaultans = {'6.25','10','80','sqrt(4/3)'};          % default inputs
% options.Resize = 'on';                                                  % allows for resizing of box
% answer = inputdlg(prompt, dlg_title, [1 50], defaultans, options);      % creates box
% epsilon = str2double(answer{1});                
% minPts = str2double(answer{2});                 
% percent = str2num(answer{3});
% s = str2num(answer{4});
%}

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
        xlabel('DAT WAY -->')
    else
        xlabel('<-- DAT WAY')
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

%% Cluster Checks
for d = 7:7
   thisImage = ClusterInfo{d,1};
   pic = ClusterInfo{d,2}; %pic = I
   
   numClust = length(thisImage);

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
   
   %Remove bottom border
   for i = 1:numClust %Iterate through Clusters
       if thisImage(i).RemoveCluster == 0
           clustPoints = thisImage(i).ClusterIndices; %Get cluster indices
           for a = 1:length(clustPoints(:,1)) %Search through cluster indices
               if (pic((clustPoints(a,2)+5), clustPoints(a,1)) == 0) %If pixel below any cluster has intensity 0, mark cluster for removal
                   thisImage(i).RemoveCluster = 1;
               end
           end       
       end 
   end    
   fprintf('Removed clusters from bottom border\n');
   
   %Rudimentary check for vessels(just based on thinness, not the real
   %vesselness check) 
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

            if(ylength > 20 || xlength > 20 || DiagnolLength > 20)
                thisImage(t).RemoveCluster = 1;
            end

            if(ylength >= xlength*3 || xlength>= ylength*3)
                thisImage(t).RemoveCluster = 1;
            end 
        end 
    end
    fprintf('Removed vessels by shape\n');
    
    ClusterInfo{d,1} = thisImage; %Save Info to ClusterInfo
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
    title(sprintf('%s - Post Spatial Analysis',ptID));
%     plot(xunit, yunit);
    plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    
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

lowchange = lowsquarechange(topx);
lowsquarechange1 = lowavesquarechange(topx);

 %% Corresponding Nipple check: Get coordinates of nipples
% You have to track the nipples over time because of the patient's shift.
% If you want to be able to compare the exact same spot over time you need
% to do this because of patient movement. Go through and click on the
% right/left nipple in each image (right first) 
figure('Name','Select Right nipple'), 
 for i = 1:15
    I1 = I_mat{i}(find(I_mat{i}>0));
    imshow(I_mat{i}, [min(I1) max(I1)])
    hold on
    xlabel('<--')
    [X_RightNip{i},Y_RightNip{i}] = ginput(1);
 end
close

questdlg('Switch sides','Switch sides','Ok','Sure','Ok')

figure('Name','Select Left Nipple'), 
for i =1:15
    I1 = I_mat{i}(find(I_mat{i}>0));
    imshow(I_mat{i}, [min(I1) max(I1)]) % gets coordinates of nipple
    hold on,               
    xlabel('-->') 
    [X_LeftNip{i},Y_LeftNip{i}] = ginput(1);
end
close

%Movement of nipple at each minute relative to n=7
%Should be 0 for all at i = 7
for i = 1:15
    XLeftchange{i} = X_LeftNip{i} - X_LeftNip{7};
    YLeftchange{i} = Y_LeftNip{i} - Y_LeftNip{7};
    XRightchange{i} = X_RightNip{i} - X_RightNip{7};
    YRightchange{i} = Y_RightNip{i} - Y_RightNip{7};
end

figure('Name','Select Middle'), imshow(I_mat{7},[min(pic_adj) max(pic_adj)]), xlabel('Select Middle')
[middle,~] = ginput(1);
close

%% Look at clusters over time
% Remove clusters that change a lot over time 
for o = 7:7
    thisImage = ClusterInfo{o,1}; %Get Current Image Info
    numClusters = length(thisImage);
%     clusterX = cell(15,numClusters);
%     clusterY = cell(15,numClusters);
    for p = 1:numClusters
        if thisImage(p).RemoveCluster == 0
            clusterIndices = thisImage(p).ClusterIndices;
%{
%             clusterY{7,p} = clusterIndices(:,2);
%             clusterX{7,p} = clusterIndices(:,1);
%             if clusterX{7,p}(1) <= middle
%                 for i = 1:15
%                     clusterX{i,p} = floor(clusterIndices(:,1) + XRightchange{i});
%                     clusterY{i,p} = floor(clusterIndices(:,2) + YRightchange{i});
%                 end     
%             else 
%                 for i = 1:15
%                     clusterX{i,p} = floor(clusterIndices(:,1) + XLeftchange{i});
%                     clusterY{i,p} = floor(clusterIndices(:,2) + YLeftchange{i});
%                 end 
%             end 
%             for t = 1:15
%                 for j = 1:length(clusterX{t,p})
%                     if clusterX{t,p}(j) < 1
%                         clusterX{t,p}(j) = 1;
%                     end 
%                     if clusterX{t,p}(j) > xmax
%                         clusterX{t,p}(j) = xmax;
%                     end 
%                     if clusterY{t,p}(j) < 1
%                         clusterY{t,p}(j) = 1;
%                     end
%                     if clusterY{t,p}(j) > ymax
%                         clusterY{t,p}(j) = ymax;
%                     end 
%                 end 
%             end 
%}
            for q = 1:15
                I1 = I_mat{q};
%{
%                 for j =1:length(clusterY{q,p})
%                     I2(j) = I1(clusterY{q,p}(j),clusterX{q,p}(j));
%                 end
%                 I3 = I2(find(I2>0));
                %averages of same cluster over time
%                 avgs(q,p) = mean2(I3);
%                 clear I2 I3;
%}
                avgs(q,p) = mean2(I1(clusterIndices(:,2), clusterIndices(:,1)));
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
            if(abs(avgStepChange(l)) > abs(lowsquarechange1)) %If the average change is too high
                thisImage(l).RemoveCluster = 1;
            end
            if(abs(totalChange(l)) < 75 && totalChange(l) > 0)
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
    title(sprintf('%s - Post Temporal Analysis',ptID));
%     plot(xunit, yunit);
    plot([c1(1) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end
%}

%% Vessel Check
% Checking for vessels. Will need to be updated to include fibermetric in
% the future. Doesn't work very well right now. 
for o = 7:7 
    thisImage = ClusterInfo{o,1}; %Get Current Image Info

    numClust = length(thisImage);
    
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
        
        if counter > 1
            PerpCentXWidth(t) = max(cell2mat(XDistances));
            PerpCentYWidth(t) = max(cell2mat(YDistances));
        else 
            PerpCentXWidth(t) = 0;
            PerpCentYWidth(t) = 0;
        end 
        
        CenX = ClusterInfo{7,1}(t).ClusterCentroid(1);
        CenY = ClusterInfo{7,1}(t).ClusterCentroid(2);
        
        counter = 1;
        for i = 1:length(thisImage(t).ClusterIndices)
            if ClusterInfo{7,1}(t).ClusterIndices(i,1) == CenX
                counter = counter + 1;
            end
        end
        CentXWidth(t) = counter;
        
        counter = 1;
        for i = 1:length(thisImage(t).ClusterIndices)
            if ClusterInfo{7,1}(t).ClusterIndices(i,2) == CenY
                counter = counter + 1;
            end
        end
        CentYWidth(t) = counter;

        
        for j = 1:length(ClusterInfo{7,1}(t).ClusterIndices(:,1))
            CentMovementXPos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) + (CentXWidth(t));
            CentMovementYPos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) + (CentYWidth(t)); 
            PerpMovementXPos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) + (PerpCentXWidth(t));
            PerpMovementYPos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) + (PerpCentYWidth(t));
            if  CentMovementXPos{j} > c
                CentMovementXPos{j} = c;
            end
            if  PerpMovementXPos{j} > c
                PerpMovementXPos{j} = c;
            end
            if  CentMovementYPos{j} > r
                CentMovementYPos{j} = r;
            end
            if  PerpMovementYPos{j} > r
                PerpMovementYPos{j} = r;
            end
        end
        
        for j = 1:length(ClusterInfo{7,1}(t).ClusterIndices(:,1))
            CentMovementXNeg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) - (CentXWidth(t));
            CentMovementYNeg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) - (CentYWidth(t));
            PerpMovementXNeg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) - (PerpCentXWidth(t));
            PerpMovementYNeg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) - (PerpCentYWidth(t));
            if CentMovementXNeg{j} < 1
                CentMovementXNeg{j} = 1;
            end
            if PerpMovementXNeg{j} < 1
                PerpMovementXNeg{j} = 1;
            end
            if CentMovementYNeg{j} < 1
                CentMovementYNeg{j} = 1;
            end
            if PerpMovementYNeg{j} < 1
                PerpMovementYNeg{j} = 1;
            end
        end
        
        for k = 1:length(ClusterInfo{7,1}(t).ClusterIndices(:,1))
            OriginalClusterYIndices{k} = ClusterInfo{7,1}(t).ClusterIndices(k,2);
            OriginalClusterXIndices{k} = ClusterInfo{7,1}(t).ClusterIndices(k,1);
        end


        NewVesClustIndicesPos{t} = transpose([CentMovementXPos;CentMovementYPos]);
        NewVesClustIndicesNeg{t} = transpose([CentMovementXNeg;CentMovementYNeg]);
        NewVesClustIndicesRight{t} = transpose([CentMovementXPos;OriginalClusterYIndices]);
        NewVesClustIndicesLeft{t} = transpose([CentMovementXNeg;OriginalClusterYIndices]);
        NewVesClustIndicesUp{t} = transpose([OriginalClusterXIndices;CentMovementYPos]);
        NewVesClustIndicesDown{t} = transpose([OriginalClusterXIndices;CentMovementYNeg]);
        OldNewVesClustIndicesPos{t} = transpose([PerpMovementXPos; PerpMovementYPos]);
        OldNewVesClustIndicesNeg{t} = transpose([PerpMovementXNeg; PerpMovementYNeg]);
        
        
        clear CentMovementXPos PerpMovementXPos CentMovementYPos PerpMovementYPos ...
            CentMovementXNeg PerpMovementXNeg CentMovementYNeg PerpMovementYNeg ...
            OriginalClusterYIndices OriginalClusterXIndices
        
        %New clusters from movement
        AdjustedVesselsPos = I_mat{7}(cell2mat(NewVesClustIndicesPos{t}(:,2)), cell2mat(NewVesClustIndicesPos{t}(:,1)));
        AdjustedVesselsNeg = I_mat{7}(cell2mat(NewVesClustIndicesNeg{t}(:,2)), cell2mat(NewVesClustIndicesNeg{t}(:,1)));
        AdjustedVesselsRight = I_mat{7}(cell2mat(NewVesClustIndicesRight{t}(:,2)), cell2mat(NewVesClustIndicesRight{t}(:,1)));
        AdjustedVesselsLeft = I_mat{7}(cell2mat(NewVesClustIndicesLeft{t}(:,2)), cell2mat(NewVesClustIndicesLeft{t}(:,1)));
        AdjustedVesselsUp = I_mat{7}(cell2mat(NewVesClustIndicesUp{t}(:,2)), cell2mat(NewVesClustIndicesUp{t}(:,1)));
        AdjustedVesselsDown = I_mat{7}(cell2mat(NewVesClustIndicesDown{t}(:,2)), cell2mat(NewVesClustIndicesDown{t}(:,1)));
        OldAdjustedVesselsPos = I_mat{7}(cell2mat(OldNewVesClustIndicesPos{t}(:,2)), cell2mat(OldNewVesClustIndicesPos{t}(:,1)));
        OldAdjustedVesselsNeg = I_mat{7}(cell2mat(OldNewVesClustIndicesNeg{t}(:,2)), cell2mat(OldNewVesClustIndicesNeg{t}(:,1)));

        %Adjust clusters to take nonzero values
        AdjVesPos = AdjustedVesselsPos(find(AdjustedVesselsPos>0));
        AdjVesNeg = AdjustedVesselsNeg(find(AdjustedVesselsNeg>0));
        AdjVesRight = AdjustedVesselsRight(find(AdjustedVesselsRight>0));
        AdjVesLeft = AdjustedVesselsLeft(find(AdjustedVesselsLeft>0));
        AdjVesUp = AdjustedVesselsUp(find(AdjustedVesselsUp>0));
        AdjVesDown = AdjustedVesselsDown(find(AdjustedVesselsDown>0));
        OldAdjVesPos = OldAdjustedVesselsPos(find(OldAdjustedVesselsPos>0));
        OldAdjVesNeg = OldAdjustedVesselsNeg(find(OldAdjustedVesselsNeg>0));

        %Find average intensity of each cluster
        avgAdjustedVesselPos(t) = mean2(AdjVesPos);
        avgAdjustedVesselNeg(t) = mean2(AdjVesNeg);
        avgAdjustedVesselRight(t) = mean2(AdjVesRight);
        avgAdjustedVesselLeft(t) = mean2(AdjVesLeft);
        avgAdjustedVesselUp(t) = mean2(AdjVesUp);
        avgAdjustedVesselDown(t) = mean2(AdjVesDown);
        OldAvgAdjustedVesselPos(t) = mean2(OldAdjVesPos);
        OldAvgAdjustedVesselNeg(t) = mean2(OldAdjVesNeg);
        
        %Check standard deviation of new clusters: 
        %If the cluster moved onto the nipple the standard deviationw would
        %be really high
        stdAdjustedVesselPos(t) = std2(AdjVesPos);
        stdAdjustedVesselNeg(t) = std2(AdjVesNeg);
        stdAdjustedVesselRight(t) = std2(AdjVesRight);
        stdAdjustedVesselLeft(t) = std2(AdjVesLeft);
        OldStdAdjustedVesselPos(t) = std2(OldAdjVesPos);
        OldStdAdjustedVesselNeg(t) = std2(OldAdjVesNeg);
        stdAdjustedVesselUp(t) = std2(AdjVesUp);
        stdAdjustedVesselDown(t) = std2(AdjVesDown);
        

        VesselTreshold = avgs(7,t)*0.015
%         if(avgAdjustedVesselPos(t) + VesselTreshold < avgs(7,t) && stdAdjustedVesselPos(t) < 100)
%             thisImage(t).RemoveCluster = 1;
%         end
%         if (avgAdjustedVesselNeg(t)+ VesselTreshold < avgs(7,t) && stdAdjustedVesselNeg(t) < 100)
%             thisImage(t).RemoveCluster = 1;
%         end        
        if ((OldAvgAdjustedVesselPos(t) + VesselTreshold < avgs(7,t)) && OldStdAdjustedVesselPos(t) < 100)
            thisImage(t).RemoveCluster = 1;
        end
        if (OldAvgAdjustedVesselNeg(t) + VesselTreshold < avgs(7,t) && OldStdAdjustedVesselNeg(t) < 100)
            thisImage(t).RemoveCluster = 1;
        end 
%         if(avgAdjustedVesselRight(t) + VesselTreshold < avgs(7,t) && stdAdjustedVesselRight(t) < 100)
%             thisImage(t).RemoveCluster = 1;
%         end 
%         if (avgAdjustedVesselLeft(t) + VesselTreshold < avgs(7,t) && stdAdjustedVesselLeft(t) < 100)
%             thisImage(t).RemoveCluster = 1;
%         end 
%         if(avgAdjustedVesselUp(t) + VesselTreshold < avgs(7,t) && stdAdjustedVesselUp(t) < 100)
%             thisImage(t).RemoveCluster = 1;
%         end 
%         if (avgAdjustedVesselDown(t) + VesselTreshold < avgs(7,t) && stdAdjustedVesselDown(t) < 100)
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
    plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end


%% Compare each cluster to the same region on the opposite breast
% If the cluster is warmer than the same region on the opposite breast keep
% it, otherwise mark for removal 
for o = 7:7
    thisImage = ClusterInfo{o,1};
    numClust = length(thisImage);
    
    for i = 1:numClust
        if thisImage(i).RemoveCluster == 0
            indices = thisImage(i).ClusterIndices;
            clusterXIndices{i} = indices(:,1);
            clusterYIndices{i} = indices(:,2);
            numPoints = length(clusterXIndices{i});
            
            %Track movement of cluster 
            %Also find distance from cluster points to reference nipple
            if clusterXIndices{i}(1) <= middle
                for j = 1:15
                    shiftedXIndices{j,i} = clusterXIndices{i} + cell2mat(XRightchange(j));
                    shiftedYIndices{j,i} = clusterYIndices{i} + cell2mat(YRightchange(j));
                    nip2ClustX{j,i} = X_RightNip{j} - shiftedXIndices{j,i};
                    nip2ClustY{j,i} = Y_RightNip{j} - shiftedYIndices{j,i};
                end
            else
                for j = 1:15
                    shiftedXIndices{j,i} = clusterXIndices{i} + cell2mat(XLeftchange(j));
                    shiftedYIndices{j,i} = clusterYIndices{i} + cell2mat(YLeftchange(j));
                    nip2ClustX{j,i} = X_RightNip{j} - shiftedXIndices{j,i};
                    nip2ClustY{j,i} = Y_RightNip{j} - shiftedYIndices{j,i};
                end
            end
            
            %Map to corresponding region
            if clusterXIndices{i}(1) > middle %if cluster was on left; cooresponding region on right
                for j = 1:15
                    corrXIndices{j,i} = X_RightNip{j} + nip2ClustX{j,i};
                    corrYIndices{j,i} = Y_RightNip{j} + nip2ClustX{j,i};
                end
            else 
                for j = 1:15
                    corrXIndices{j,i} = X_LeftNip{j} + nip2ClustX{j,i};
                    corrYIndices{j,i} = Y_LeftNip{j} + nip2ClustX{j,i};
                end
            end
            
            
            for j = 1:15
                for x = 1:numPoints
                    if shiftedXIndices{j,i}(x) < 1
                        shiftedXIndices{j,i}(x) = 1;
                    end
                    if shiftedXIndices{j,i}(x) > c 
                        shiftedXIndices{j,i}(x) = c;
                    end 
                    if shiftedYIndices{j,i}(x) < 1
                        shiftedYIndices{j,i}(x) = 1;
                    end
                    if shiftedYIndices{j,i}(x) > r 
                        shiftedYIndices{j,i}(x) = r;
                    end
                    if corrXIndices{j,i}(x) < 1
                        corrXIndices{j,i}(x) = 1;
                    end
                    if corrXIndices{j,i}(x) > c
                        corrXIndices{j,i}(x) = c;
                    end 
                    if corrYIndices{j,i}(x) < 1
                        corrYIndices{j,i}(x) = 1;
                    end
                    if corrYIndices{j,i}(x) > r 
                        corrYIndices{j,i}(x) = r;
                    end
                end
                clustIntensities(x) = I_mat{j}(floor(shiftedYIndices{j,i}(x)),floor(shiftedXIndices{j,i}(x)));
                corrIntensities(x)= I_mat{j}(floor(corrYIndices{j,i}(x)),floor(corrXIndices{j,i}(x)));
                
                clustAverageTracked{i,j} = mean2(clustIntensities);
                corrAverageTracked{i,j} = mean2(corrIntensities);
            end
            
            clusterChangeTotal(i) = clustAverageTracked{i,15} - clustAverageTracked{i,1};
            corrChangeTotal(i) = corrAverageTracked{i,15} - corrAverageTracked{i,1};
            for y = 1:14
                clusterChangeStep(i,y) = clustAverageTracked{i,y+1} - clustAverageTracked{i,y};
                corrChangeStep(i,y) = corrAverageTracked{i,y+1} - corrAverageTracked{i,y};
            end 
            clusterAvgStepChange(i) = mean2(clusterChangeStep(i,:));
            corrAvgStepChange(i) = mean2(corrChangeStep(i,:));
            
            if abs(clusterChangeTotal(i)) > abs(corrChangeTotal(i))+75 && clusterChangeTotal(i) < 8000
                thisImage(i).RemoveCluster = 1;
            end 
            if clusterChangeTotal(i) > 0 && clusterChangeTotal(i) < 8000
                thisImage(i).RemoveCluster = 1;
            end
            
        end
    end 
end 

ClusterInfo{o,1} = thisImage;

%% Plot results (again) 
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
    title(sprintf('%s - Post Corresponding Region Check',ptID));
%     plot(xunit, yunit);
    plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    

   
    ClusterInfo{7,1} = thisImage;
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end
fprintf('DBSCAN_3D took %04f seconds to run\n',toc)
