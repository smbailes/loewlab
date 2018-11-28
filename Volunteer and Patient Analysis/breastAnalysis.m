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
    
figure('Name','Select nipple (w/o Tumor)'), 
for i = 1:15 % Get coordinates of both nipples at each minute
    I1 = I_mat{i};              % Display first image
    I = getMatrixOutliers(I1);  % Remove outliers
    I_adj1 = I(find(I>0));
    imshow(I1, [min(I_adj1) max(I_adj1)]) % gets coordinates of nipple
    hold on,
    if strcmp(sideString,'Left') == 1                   % Direct user to tumor side
        xlabel('<--')
    else
        xlabel('-->')
    end 
    [XCorrNip{i},YCorrNip{i}] = ginput(1);
end

questdlg('Switch sides','Switch sides','Ok','Sure','Ok')

figure('Name','Select Nipple (w/ Tumor)'), 
for i =1:15
    I1 = I_mat{i};              % Display first image
    I = getMatrixOutliers(I1);  % Remove outliers
    I_adj1 = I(find(I>0));
    imshow(I1, [min(I_adj1) max(I_adj1)]) % gets coordinates of nipple
    hold on,
    if strcmp(sideString,'Left') == 1                   % Direct user to tumor side
        xlabel('-->')
    else
        xlabel('<--')
    end 
    [XTumNip{i},YTumNip{i}] = ginput(1); 
end 
     
% Get a center of tumor for each time
for i = 1:15
    XTumNew{i} = XTumNip{i} + dx;                       % Coordinates of center of tumor
    YTumNew{i} = YTumNip{i} - dy; 
end
    xbox = xbox * scale;                    % Convert X and Y dimensions of tumor from cm to pixels
    ybox = ybox * scale;
    xbox = xbox*1.2;
    ybox = ybox*1.2;
    
    cT1 = [round(XTumNew{7} - xbox/2), round(YTumNew{7} - ybox/2)];  % Top Left Corner
    cT2 = [round(XTumNew{7} - xbox/2), round(YTumNew{7} + ybox/2)];  % Bottom Left Corner
    cT3 = [round(XTumNew{7} + xbox/2), round(YTumNew{7} + ybox/2)];  % Bottum Right Corner
    cT4 = [round(XTumNew{7} + xbox/2), round(YTumNew{7} - ybox/2)];  % Top Right Corner
%     th = 0:pi/50:2*pi;
%     if xbox <= ybox
%         xunit = xbox * cos(th) + XTumNew{7};
%         yunit = xbox * sin(th) + YTumNew{7};
%     elseif ybox < xbox
%         xunit = ybox * cos(th) + XTumNew{7};
%         yunit = ybox * sin(th) + YTumNew{7};      
%     end 
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
    c1 = [round(XCorrNew{7} - xbox/2), round(YCorrNew{7} - ybox/2)];  % Top Left Corner
    c2 = [round(XCorrNew{7} - xbox/2), round(YCorrNew{7} + ybox/2)];  % Bottom Left Corner
    c3 = [round(XCorrNew{7} + xbox/2), round(YCorrNew{7} + ybox/2)];  % Bottum Right Corner
    c4 = [round(XCorrNew{7} + xbox/2), round(YCorrNew{7} - ybox/2)];  % Top Right Corner      
%     th = 0:pi/50:2*pi;
%     if xbox <= ybox
%         xunit1 = xbox * cos(th) + XCorrNew{7};
%         yunit1 = xbox * sin(th) + YCorrNew{7};
%     elseif ybox < xbox
%         xunit1 = ybox * cos(th) + XCorrNew{7};
%         yunit1 = ybox * sin(th) + YCorrNew{7};      
%     end 
    hold off  
    close    

%% Select ROI - tumor
figure('Name', 'Select Tumor Region'),
imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast
hold on;
plot([cT1(1) cT2(1)],[cT1(2) cT2(2)],'b');                      % Create red box region on Image Display
plot([cT2(1) cT3(1)],[cT2(2) cT3(2)],'b');
plot([cT3(1) cT4(1)],[cT3(2) cT4(2)],'b');
plot([cT4(1) cT1(1)],[cT4(2) cT1(2)],'b');


e = imrect();

xyTum = wait(e);
hold on;
%% Select ROI - corresponding region
figure('Name', 'Select Coresponding Region'),
imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast
hold on;
plot([c1(1 ) c2(1)],[c1(2) c2(2)],'b');                      % Create red box region on Image Display
plot([c2(1) c3(1)],[c2(2) c3(2)],'b');
plot([c3(1) c4(1)],[c3(2) c4(2)],'b');
plot([c4(1) c1(1)],[c4(2) c1(2)],'b');
e = imrect();

xyCorr = wait(e);
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

for i = 1:15
    TumorStdev{i} = std2((TumRegion1{i}));
    CorrStdev{i} = std2((CorrRegion1{i}));
end 

for i = 1:14
    TumStepChange{i} = TumAve{i+1} - TumAve{i};
    CorrStepChange{i} = CorrAve{i+1} - CorrAve{i};
end

tumAveInitial = TumAve{1}, tumAveFinal = TumAve{15}, 
corrAveInitial = CorrAve{1}, corrAveFinal = CorrAve{15}
tumStdInitial = TumorStdev{1}, tumStdFinal = TumorStdev{15},
corrStdIntial = CorrStdev{1}, corrStdFinal = CorrStdev{15},

TumAveStepChange = mean2(cell2mat(TumStepChange))
CorrAveStepChange = mean2(cell2mat(CorrStepChange))
TumTotChange = TumAve{15} - TumAve{1}
CorrTotChange = CorrAve{15} - CorrAve{1}


for i = 1:15
    if TumAve{i} < CorrAve{i}
        fprintf('Corresponding region is warmer at n=%f \n', i)
    end
end

close all
%% Plot them together
figure, 
% subplot(2,1,1),
% title('Tumor Region')
imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast
hold on;
p1 = plot([cT1(1) cT2(1)],[cT1(2) cT2(2)],'b');                      % Create red box region on Image Display
p2 = plot([cT2(1) cT3(1)],[cT2(2) cT3(2)],'b');
p3 = plot([cT3(1) cT4(1)],[cT3(2) cT4(2)],'b');
p4 = plot([cT4(1) cT1(1)],[cT4(2) cT1(2)],'b');

subplot(2,1,2),
title('Corresponding Region')
imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast
hold on;
p5 = plot([c1(1 ) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
p6 = plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
p7 = plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
p8 = plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
    


%% Plot temperature changes
    
for i = 1:15
    tumorAverage(i) = TumAve{i};
    corrAverage(i) = CorrAve{i};
end 

figure, plot(1:15, tumorAverage,'o-')
hold on
plot(1:15, corrAverage,'o-')
legend('Tumor Region','Corresponding Region')
xlabel('Time (minutes)')
ylabel('Temperature (Digital Counts)')