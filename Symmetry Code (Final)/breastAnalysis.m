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
    I_sort1 = sort(I_adj1);
    figure('Name','Nipple Identification (side with tumor)')
    imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast

    if strcmp(sideString,'Left') == 1                   % Direct user to tumor side
        xlabel('-->')
    else
        xlabel('<--')
    end 
    hold on
    fprintf('Select Reference Nipple (side with tumor)\n');  % User input of nipple region
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

    figure('Name','Nipple Identification (side without tumor)')
    imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast

    if strcmp(sideString,'Left') == 1                   % Direct user to tumor side
        xlabel('<--')
    else
        xlabel('-->')
    end 
    hold on
    fprintf('Select Reference Nipple (side without tumor)\n');  % User input of nipple region
    [Xorg1,Yorg1] = ginput(1);
    plot(Xorg1,Yorg1, '*');                   % Plot center of nipple
    Xnew1 = Xorg1 + dx1;                       % Coordinates of center of tumor
    Ynew1 = Yorg1 - dy1;    
 
    th = 0:pi/50:2*pi;
    if xbox <= ybox
        xunit1 = xbox * cos(th) + Xnew1;
        yunit1 = xbox * sin(th) + Ynew1;
    elseif ybox < xbox
        xunit1 = ybox * cos(th) + Xnew1;
        yunit1 = ybox * sin(th) + Ynew1;      
    end 
    hold off  
    close    

%% Select ROI - tumor
figure('Name', 'Select Tumor Region'),
imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast
hold on;
plot(xunit, yunit);
e = imellipse();

xy = wait(e);
binaryImage = e.createMask();
BW_t = uint16(binaryImage);
hold on;

I1 = I_mat{n};
I2 = I_mat{n}(find(I_mat{n}>0));

I3 = I1.*BW_t;
tumorRegion = I3(find(I3>0));

%% Select ROI - corresponding region
figure('Name', 'Select Coresponding Region'),
imshow(I,[min(I_adj1) max(I_adj1)]);               % Display with contrast
hold on;
plot(xunit1, yunit1);
e = imellipse();

xy = wait(e);
binaryImage = e.createMask();
BW_c = uint16(binaryImage);
hold on;

I1 = I_mat{n};
I2 = I_mat{n}(find(I_mat{n}>0));

I3 = I1.*BW_c;
corrRegion = I3(find(I3>0));


%% Show histograms

figure('Name', 'Histogram');
    I1 = I_mat{n};
    I2 = I_mat{n}(find(I_mat{n}>0));

%     subplot(4,4,n)

    histogram(I2,500,'FaceColor','r','EdgeColor','r');
    title('Entire Image(red) v. Tumor Region(yellow) v. Corresponding Region(blue)')
    hold on
    yyaxis right
    ylim([0 100])
    histogram(tumorRegion,500,'FaceColor','y','EdgeColor','y');
    hold on
    histogram(corrRegion, 500, 'FaceColor', 'b', 'EdgeColor', 'b');
    
%% Find averages

for i = 1:15
    
    I = I_mat{i};
    
    I_t = I.*BW_t;
    I_c = I.*BW_c;
    
    I_tumor = I_t(find(I_t>0));
    I_corr = I_c(find(I_c>0));
    
    tumorAv(i) = mean2(I_tumor);
    corrAv(i) = mean2(I_corr);
    
    diff(i) = abs(tumorAv(i)-corrAv(i));
    
end


    

    
