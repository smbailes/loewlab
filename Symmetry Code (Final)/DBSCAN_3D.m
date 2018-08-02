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
% close all;

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
    s = sqrt(10/8);
elseif range(1) <= 2700
    s = 1;
end     

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
            
            if(ylength < 5 || xlength < 5)
                thisImage(t).RemoveCluster = 1;
            end

            if(ylength >= xlength*5 || xlength>= ylength*5)
                thisImage(t).RemoveCluster = 1;
            end 
        end 
    end
    fprintf('Removed vessels by shape\n');
    
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

 %% Corresponding Nipple check: Get coordinates of nipples

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
%Movement of nipple at each minute
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

for o = 7:7
    thisImage = ClusterInfo{o,1}; %Get Current Image Info
    numClusters = length(thisImage);
    
    for p = 1:numClusters
        if thisImage(p).RemoveCluster == 0
            clusterIndices = thisImage(p).ClusterIndices;
            clusterX = cell(15,1);
            clusterY = cell(15,1);
            clusterX{7} = clusterIndices(:,2);
            clusterY{7} = clusterIndices(:,1);
            if clusterX{7}(1) <= middle
                for i = 1:15
                    clusterX{i} = floor(clusterX{7} + XRightchange{i});
                    clusterY{i} = floor(clusterY{7} + YRightchange{i});
                end 
            else 
                for i = 1:15
                    clusterX{i} = floor(clusterX{7} + XLeftchange{i});
                    clusterY{i} = floor(clusterY{7} + YLeftchange{i});
                end 
            end 
            for q = 1:15
                I1 = I_mat{q};
                %averages of same cluster over time
                avgs(q,p) = mean2(I1(clusterX{i},clusterY{i}));
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
        
        OldXWidth(t) = max(cell2mat(XDistances));
        OldYWidth(t) = max(cell2mat(YDistances));
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
            OldNewVesClustXpos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) + (OldXWidth(t)+5);
            OldNewVesClustYpos{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) + (OldYWidth(t)+5);
            if  NewVesClustXpos{j} > Xdimen
                NewVesClustXpos{j} = Xdimen;
            end
            if  OldNewVesClustXpos{j} > Xdimen
                OldNewVesClustXpos{j} = Xdimen;
            end
            if  NewVesClustYpos{j} > Ydimen
                NewVesClustYpos{j} = Ydimen;
            end
            if  OldNewVesClustYpos{j} > Ydimen
                OldNewVesClustYpos{j} = Ydimen;
            end
        end
        
        for j = 1:length(ClusterInfo{7,1}(t).ClusterIndices(:,1))
            NewVesClustXneg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) - (XWidth(t)+5);
            NewVesClustYneg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) - (YWidth(t)+5);
            OldNewVesClustXneg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,1) - (OldXWidth(t)+5);
            OldNewVesClustYneg{j} = ClusterInfo{7,1}(t).ClusterIndices(j,2) - (OldYWidth(t)+5);
            if NewVesClustXneg{j} < 1
                NewVesClustXneg{j} = 1;
            end
            if OldNewVesClustXneg{j} < 1
                OldNewVesClustXneg{j} = 1;
            end
            if NewVesClustYneg{j} < 1
                NewVesClustYneg{j} = 1;
            end
            if OldNewVesClustYneg{j} < 1
                OldNewVesClustYneg{j} = 1;
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
        OldNewVesClustIndicesPos{t} = transpose([OldNewVesClustXpos; OldNewVesClustYpos]);
        OldNewVesClustIndicesNeg{t} = transpose([OldNewVesClustXneg; OldNewVesClustYneg]);
        
        
        clear NewVesClustXpos OldNewVesClustXpos NewVesClustYpos OldNewVesClustYpos ...
            NewVesClustXneg OldNewVesClustXneg NewVesClustYneg OldNewVesClustYneg ...
            OriginalClusterYIndices OriginalClusterXIndices
        
        AdjustedVesselsPos = I_mat{7}(cell2mat(NewVesClustIndicesPos{t}(:,2)), cell2mat(NewVesClustIndicesPos{t}(:,1)));
        AdjustedVesselsNeg = I_mat{7}(cell2mat(NewVesClustIndicesNeg{t}(:,2)), cell2mat(NewVesClustIndicesNeg{t}(:,1)));
        AdjustedVesselsRight = I_mat{7}(cell2mat(NewVesClustIndicesRight{t}(:,2)), cell2mat(NewVesClustIndicesRight{t}(:,1)));
        AdjustedVesselsLeft = I_mat{7}(cell2mat(NewVesClustIndicesLeft{t}(:,2)), cell2mat(NewVesClustIndicesLeft{t}(:,1)));
        OldAdjustedVesselsPos = I_mat{7}(cell2mat(OldNewVesClustIndicesPos{t}(:,2)), cell2mat(OldNewVesClustIndicesPos{t}(:,1)));
        OldAdjustedVesselsNeg = I_mat{7}(cell2mat(OldNewVesClustIndicesNeg{t}(:,2)), cell2mat(OldNewVesClustIndicesNeg{t}(:,1)));

        AdjVesPos = AdjustedVesselsPos(find(AdjustedVesselsPos>0));
        AdjVesNeg = AdjustedVesselsNeg(find(AdjustedVesselsNeg>0));
        AdjVesRight = AdjustedVesselsRight(find(AdjustedVesselsRight>0));
        AdjVesLeft = AdjustedVesselsLeft(find(AdjustedVesselsLeft>0));
        OldAdjVesPos = OldAdjustedVesselsPos(find(OldAdjustedVesselsPos>0));
        OldAdjVesNeg = OldAdjustedVesselsNeg(find(OldAdjustedVesselsNeg>0));


        avgAdjustedVesselPos(t) = mean2(AdjVesPos);
        avgAdjustedVesselNeg(t) = mean2(AdjVesNeg);
        avgAdjustedVesselRight(t) = mean2(AdjVesRight);
        avgAdjustedVesselLeft(t) = mean2(AdjVesLeft);
        OldAvgAdjustedVesselPos(t) = mean2(OldAdjVesPos);
        OldAvgAdjustedVesselNeg(t) = mean2(OldAdjVesNeg);
        
        stdAdjustedVesselPos(t) = std2(AdjVesPos);
        stdAdjustedVesselNeg(t) = std2(AdjVesNeg);
        stdAdjustedVesselRight(t) = std2(AdjVesRight);
        stdAdjustedVesselLeft(t) = std2(AdjVesLeft);
        OldStdAdjustedVesselPos(t) = std2(OldAdjVesPos);
        OldStdAdjustedVesselNeg(t) = std2(OldAdjVesNeg);
        

        VesselPDiff = avgs(7,t)*0.03
%         if(avgAdjustedVesselPos(t) + VesselPDiff < avgs(7,t) && stdAdjustedVesselPos(t) < 100)
%             thisImage(t).RemoveCluster = 1;
%         end
%         if (avgAdjustedVesselNeg(t)+ VesselPDiff < avgs(7,t) && stdAdjustedVesselNeg(t) < 100)
%             thisImage(t).RemoveCluster = 1;
%         end        
        if ((OldAvgAdjustedVesselPos(t) + VesselPDiff < avgs(7,t)) && OldStdAdjustedVesselPos(t) < 100)
            thisImage(t).RemoveCluster = 1;
        end
        if (OldAvgAdjustedVesselNeg(t) + VesselPDiff < avgs(7,t) && OldStdAdjustedVesselNeg(t) < 100)
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
    plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
    plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
    plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
    plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    hold on 
    
    ClusterInfo{g,3} = clustData; %Save updated Cluster Info to Array
end


%% Track clusters over time
% numClust = length(ClustStruct);

clear ClustInfoCell RemoveCluster NumClust JustClust Indices XClusterIndices YClusterIndices...
    XLeftchange YLeftchange XRightchange YRightchange AdjustedClustStruct Values Xnip2Clust Ynip2Clust...
    XCorrIndices YCorrIndices RemainingClust XCorrIndice YCorrIndice
ClustInfoCell = struct2cell(ClusterInfo{7,1}); %converts data into cell array
RemoveCluster = cell2mat(ClustInfoCell(7,:,:)); %separates the data indicating to remove clusters
counter = 0;
for i = 1:length(ClustInfoCell)
    if RemoveCluster(i) == 0
        counter = counter+1;
        JustClust{counter} = ClustInfoCell(:,:,i);
    end
end
   RemainingClust = numel(JustClust); %finds number of clusters
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
    XLeftchange{i} = X_LeftNip{i} - X_LeftNip{7};
    YLeftchange{i} = Y_LeftNip{i} - Y_LeftNip{7};
    XRightchange{i} = X_RightNip{i} - X_RightNip{7};
    YRightchange{i} = Y_RightNip{i} - Y_RightNip{7};
end

AdjustedClustStruct = struct('ClusterNumber',d,'TumorBreastClustorXPoints',XClusterIndices,'TumorBreastClustorYPoints',YClusterIndices);

NewClustXPoints = cell(15,RemainingClust);
NewClustYPoints = cell(15,RemainingClust);
[totRow,totCol] = size(I_mat{7});

for i = 1:RemainingClust
    if mean2(XClusterIndices{i}) <= middle
        for j = 1:15
        NewClustXPoints{j,i} = XClusterIndices{i} + cell2mat(XRightchange(j)); % adjusts indices by the change in the nipple relative to time 7
        NewClustYPoints{j,i} = YClusterIndices{i} + cell2mat(YRightchange(j));
        end % New Xpoints. i = number cluster and j = points at time
    else
        for j = 1:15
        NewClustXPoints{j,i} = XClusterIndices{i} + cell2mat(XLeftchange(j));
        NewClustYPoints{j,i} = YClusterIndices{i} + cell2mat(YLeftchange(j));
        end
    end
end
[r,c] = size(NewClustXPoints);
for i = 1:RemainingClust
    if mean2(XClusterIndices{i}) <= middle
        for j = 1:15
         Xnip2Clust{j,i} = X_RightNip{j} - NewClustXPoints{j,i};
         Ynip2Clust{j,i} = Y_RightNip{j} - NewClustYPoints{j,i};
        end
    else
        for j = 1:15
         Xnip2Clust{j,i} = X_LeftNip{j} - NewClustXPoints{j,i};
         Ynip2Clust{j,i} = Y_LeftNip{j} - NewClustYPoints{j,i};
        end
    end
end
for i = 1:RemainingClust
    if mean2(XClusterIndices{i}) <= middle
        for j = 1:15 % 
            XCorrIndices{j,i} = X_LeftNip{j} + Xnip2Clust{j,i};
            YCorrIndices{j,i} = Y_LeftNip{j} + Ynip2Clust{j,i};
        end
    else
        for j = 1:15
            XCorrIndices{j,i} = X_RightNip{j} + Xnip2Clust{j,i};
            YCorrIndices{j,i} = Y_RightNip{j} + Ynip2Clust{j,i};
        end
    end
end
for j = 1:c % cluster
    for i = 1:r % time
    XClustIndice = cell2mat(NewClustXPoints(i,j)); %Creates the points for the cluster at this time
    YClustIndice = cell2mat(NewClustYPoints(i,j));
    XCorrIndice = cell2mat(XCorrIndices(i,j));
    YCorrIndice = cell2mat(YCorrIndices(i,j));
    L = length(XClustIndice);
      for k = 1:L % pixel
          if XClustIndice(k) <1 % adjusts indices if exceed image
              XClustIndice(k) = 1;              
          elseif XClustIndice(k) > totCol
              XClustIndice(k) = totCol;  
          end
          if XCorrIndice(k) <1
              XCorrIndice(k) = 1;
          elseif XCorrIndice(k) > totCol
              XCorrIndice(k) = totCol;
          end
          if YClustIndice(k) <1
              YClustIndice(k) = 1;
          elseif YClustIndice(k) > totRow
              YClustIndice(k) = totRow;
          end
          if YCorrIndice(k) <1
              YCorrIndice(k) = 1;
          elseif YCorrIndice(k) > totRow
              YCorrIndice(k) = totRow;
          end
          
          ClustValues{k} = I_mat{i}(floor(YClustIndice(k)),floor(XClustIndice(k))); %Records the value at each time
          CorrValues{k} = I_mat{i}(floor(YCorrIndice(k)),floor(XCorrIndice(k)));
      end
      
      ClustValues = cell2mat(ClustValues);
      CorrValues = cell2mat(CorrValues); 
      
      TimeClusterData{i,j} = mean(ClustValues(find(ClustValues>0))); % Records the average cluster value at each time
      CorrData{i,j} = mean(CorrValues(find(CorrValues>0))); %Records average cluster value for region opposite cluster at each time
      

      clear ClustValues CorrValues
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
       if RemoveClusterCounter <= length(NumberOfRemoval)
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

