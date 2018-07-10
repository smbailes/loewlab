clear all, close all
%% Patient Selection
    [location, ptID] = pathfinder; 

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end
    
    n=7; %Use the7th image for now

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

    
%% ROI Identification on First Image - Tumor side

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
    
    I1 = I_mat{7};              % Display first image
    I = getMatrixOutliers(I1);  % Remove outliers
    I_adj1 = I1(find(I1>0));    % Remove zero pixels
    I_sort1 = sort(I_adj1)
    %% Select Both Nipples at Each Minute
for i = 1:15 % Get coordinates of both nipples at each minute
    figure('Name','Select nipple (w/o Tumor)'), imshow(I_mat{i}, [min(I_adj1) max(I_adj1)]) % gets coordinates of nipple
    [XCorrNip{i},YCorrNip{i}] = ginput(1), close
end
questdlg('Switch sides','Switch sides','Ok','Sure','Ok')

for i =1:15
    figure('Name','Select Nipple (w/ Tumor)'), imshow(I_mat{i} , [min(I_adj1) max(I_adj1)]) % gets coordinates of nipple
    [XTumNip{i},YTumNip{i}] = ginput(1), close
end 
     
% Get a center of tumor for each time
for i = 1:15
    XTumNew{i} = XTumNip{i} + dx;                       % Coordinates of center of tumor
    YTumNew{i} = YTumNip{i} - dy; 
end
    xbox = xbox * scale;                    % Convert X and Y dimensions of tumor from cm to pixels
    ybox = ybox * scale;    
    th = 0:pi/50:2*pi;
    if xbox <= ybox
        xunit = xbox * cos(th) + XTumNew{7};
        yunit = xbox * sin(th) + YTumNew{7};
    elseif ybox < xbox
        xunit = ybox * cos(th) + XTumNew{7};
        yunit = ybox * sin(th) + YTumNew{7};      
    end 
    hold off  
    close    
%% ROI Identification on First Image - Corresponding side
    
    if hr==12
        hr1 = 12;
    else
        hr1 = 12-hr;
    end
    
    if hr1 <= 9 && hr1 > 3            
        hr_ang1 = (abs((hr1 - 9)) * (pi/6)) + pi;
    elseif hr1 <= 3 
        hr_ang1 = (3 - hr1) * (pi/6);
    elseif hr1 > 9
        hr_ang1 = pi - ((hr1 - 9) * (pi/6));
    end
    theta1 = hr_ang1;     % angle in radians
    
    [dx1,dy1] = pol2cart(theta1, rho); % Convert tumor location as angle & dist to pixel location 
    for i = 1:15
    XCorrNew{i} = XCorrNip{i} + dx1;
    YCorrNew{i} = YCorrNip{i} - dy1;  
    end
    th = 0:pi/50:2*pi;
    if xbox <= ybox
        xunit1 = xbox * cos(th) + XCorrNew{7};
        yunit1 = xbox * sin(th) + YCorrNew{7};
    elseif ybox < xbox
        xunit1 = ybox * cos(th) + XCorrNew{7};
        yunit1 = ybox * sin(th) + YCorrNew{7};      
    end 
    hold off  
    close    

%% Select ROI - tumor
figure('Name', 'Select Tumor Region'),
imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast
hold on;
plot(xunit, yunit);
e = imrect();

xyTum = wait(e);
hold on;
%% Select ROI - corresponding region
figure('Name', 'Select Coresponding Region'),
imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast
hold on;
plot(xunit1, yunit1);
e = imrect();

xyCorr = wait(e);
binaryImage = e.createMask();
BW_c = uint16(binaryImage);
hold on;
%% Show histograms

% figure('Name', 'Histogram');
%     I1 = I_mat{n};
%     I2 = I_mat{n}(find(I_mat{n}>0));
% 
% %     subplot(4,4,n)
% 
%     histogram(I2,500,'FaceColor','r','EdgeColor','r');
%     title('Entire Image(red) v. Tumor Region(yellow) v. Corresponding Region(blue)')
%     hold on
%     yyaxis right
%     ylim([0 100])
%     histogram(tumorRegion,500,'FaceColor','y','EdgeColor','y');
%     hold on
%     histogram(corrRegion, 500, 'FaceColor', 'b', 'EdgeColor', 'b');
%     
%% Find averages
for i = 1:15
    XTumChange{i} = XTumNip{7} - XTumNip{1};
    YTumChange{i} = YTumNip{7} - YTumNip{1};
    XCorrChange{i} = XCorrNip{7} - XCorrNip{1};
    YCorrChange{i} = YCorrNip{7} - YCorrNip{1};
    XTumRectInitial{i} = xyTum(1)- XTumChange{i};
    YTumRectInitial{i} = xyTum(2) - YTumChange{i};
    XCorrRectInitial{i} = xyCorr(1) - XCorrChange{i};
    YCorrRectInital{i} = xyCorr(2) - YCorrChange{i};
end
TumWidth = xyTum(3);
TumHeight = xyTum(4);
CorrWidth = xyCorr(3);
CorrHeight = xyCorr(4);

for i = 1:15
    TumRegion{i} = imcrop(I_mat{i}, [XTumRectInitial{i},YTumRectInitial{i},TumWidth,TumHeight]);
    CorrRegion{i} = imcrop(I_mat{i},[XCorrRectInitial{i},YCorrRectInital{i},CorrWidth,CorrHeight]);
    
    TumRegion1{i} = TumRegion{i}(find(TumRegion{i}>0));
    CorrRegion1{i} = CorrRegion{i}(find(CorrRegion{i}>0));
    
end
for i = 1:15
    TumAve{i} = mean2((TumRegion1{i}));
    CorrAve{i} = mean2((CorrRegion1{i}));
end 

for i = 1:14
    TumStepChange{i} = TumAve{i+1} - TumAve{i};
    CorrStepChange{i} = CorrAve{i+1} - CorrAve{i};
end
TumAve{1}, TumAve{15}, CorrAve{1}, CorrAve{15}
TumAveStepChange = mean2(cell2mat(TumStepChange))
CorrAveStepChange = mean2(cell2mat(CorrStepChange))
TumTotChange = TumAve{15} - TumAve{1}
CorrTotChange = CorrAve{15} - CorrAve{1}



    

    
