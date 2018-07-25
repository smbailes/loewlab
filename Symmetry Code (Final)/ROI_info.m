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
    notes = txt(index,7);                   % Get any notes from txt file
    celldisp(notes); 
    
    
    n = 2;
    
%% ROI Identification
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
    
    
    I_ref = I_mat{n};              % Display first image
    I = getMatrixOutliers(I_ref);  % Remove outliers
    I_adj = I_ref(find(I_ref>0));    % Remove zero pixels
    I_sort1 = sort(I_adj);
    
    [r c] = size(I_ref);
    figure('Name','ROI Identification')
    imshow(I,[min(I_adj) max(I_adj)]);               % Display with contrast
    
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
    
    th = 0:pi/50:2*pi;
    if xbox <= ybox
        xunit = xbox * cos(th) + Xnew;
        yunit = xbox * sin(th) + Ynew;
    elseif ybox < xbox
        xunit = ybox * cos(th) + Xnew;
        yunit = ybox * sin(th) + Ynew;
    end 
    plot(xunit, yunit);
    hold off  
    close    

%% Crop Circilar area
%     imshow(I,[min(I_adj) max(I_adj)]);               % Display with contrast
%     hold on;
%     plot(xunit, yunit);
%     e = imellipse();
%     xy = wait(e);
%     binaryImage = e.createMask();
%     BW = uint16(binaryImage);
%% Show Histogram 
figure('Name','Histograms over time');
subplot(5,3,1);
for n = 1:15
    I2 = I_mat{n}(find(I_mat{n}>0));
    subplot(5,3,n)
    histogram(I2,500,'FaceColor','r','EdgeColor','r');
end
%% Show circular area ROI
again = 'Yes';
figure('Name', 'Select ROI'),
imshow(I,[min(I_adj) max(I_adj)]);               % Display with contrast
hold on;
plot([c1(1 ) c2(1)],[c1(2) c2(2)],'b');                      % Create red box region on Image Display
plot([c2(1) c3(1)],[c2(2) c3(2)],'b');
plot([c3(1) c4(1)],[c3(2) c4(2)],'b');
plot([c4(1) c1(1)],[c4(2) c1(2)],'b');
e = imellipse();

while strcmp(again, 'Yes') == 1
        xy = wait(e);
        binaryImage = e.createMask();
        BW = uint16(binaryImage);
        hold on;

        figure('Name','Histogram (with ROI highlighted)');
        % subplot(4,4,1);
%         for n = 1:1

            I1 = I_mat{n};
            I2 = I_mat{n}(find(I_mat{n}>0));

            I3 = I1.*BW;
            I4 = I3(find(I3>0));

        %     subplot(4,4,n)

            histogram(I2,500,'FaceColor','r','EdgeColor','r');
            hold on
%             yyaxis right
            histogram(I4,500,'FaceColor','k','EdgeColor','k');
%         end
        
        again = questdlg('Try again?', 'Yes', 'No');

end

%% Crop Rectangular area
%     imshow(I,[min(I_adj) max(I_adj)]);               % Display with contrast
%     hold on;
%     plot(xunit, yunit);
% %     plot([c1(1) c2(1)],[c1(2) c2(2)],'r');                      % Create red box region on Image Display
% %     plot([c2(1) c3(1)],[c2(2) c3(2)],'r');
% %     plot([c3(1) c4(1)],[c3(2) c4(2)],'r');
% %     plot([c4(1) c1(1)],[c4(2) c1(2)],'r');
%     hold on;
%     hFH = imellipse();
%     xy = wait(hFH);
%     binaryImage = hFH.createMask();
%     xy = hFH.getPosition;
%     close
%% Show Rectangular Image and ROI
% figure('Name','Histograms (with ROI highlighted)');
% subplot(4,4,1);
% for n = 1:14
%     I2 = I_mat{n}(find(I_mat{n}>0));
%     newCrop = imcrop(I_mat{n}, xy);
% 
%     subplot(4,4,n)
%     histogram(I2,1000,'FaceColor','r','EdgeColor','r');
%     hold on
%     histogram(newCrop,1000,'FaceColor','k','EdgeColor','k');
% end

