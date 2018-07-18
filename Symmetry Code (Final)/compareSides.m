%Compare two sides
%after running it once for a patient comment out Location and Select
%nipples sections

clearvars -except I_mat xR yR xL yL,
close all, 
%% Location input and Nipple selection
%Only redo if the patient/volunteer has changed
redo = questdlg('New Patient/Volunteer?', 'Yes', 'No');
if strcmp(redo, 'Yes') == 1%% Input
    path  = uigetdir;
    location = strcat(path, '\');
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end
    figure('Name','Select nipple (Right)'), 
     for i = 1:15
        I1 = I_mat{i}(find(I_mat{i}>0));
        imshow(I_mat{i}, [min(I1) max(I1)])
        hold on
        xlabel('<--')
        [xR{i},yR{i}] = ginput(1);
     end
    close
    figure('Name','Select nipple (Left)'), 
    for i =1:15
        I1 = I_mat{i}(find(I_mat{i}>0));
        imshow(I_mat{i}, [min(I1) max(I1)]) % gets coordinates of nipple
        hold on,
        xlabel('-->')
        [xL{i},yL{i}] = ginput(1);
    end
    close
end

%% Parameters

prompt = {'Clock Hour', 'Distance', 'Square size'};  
dlg_title = 'Select a Region';                                         % box title
defaultans = {'12','2','2'};          % default inputs
options.Resize = 'on';                                                  % allows for resizing of box
answer = inputdlg(prompt, dlg_title, [1 50], defaultans, options);      % creates box

hrRight = str2num(answer{1});
dist = str2num(answer{2});
sqSize = str2num(answer{3});

%% Create boxes
if hrRight <= 9 && hrRight > 3            
    hr_ang = (abs((hrRight - 9)) * (pi/6)) + pi;
elseif hrRight <= 3 
    hr_ang = (3 - hrRight) * (pi/6);
elseif hrRight > 9
    hr_ang = pi - ((hrRight - 9) * (pi/6));
end
thetaRight = hr_ang;     % angle in radians

if hrRight == 12
    hrLeft = 12;
else
    hrLeft = 12-hrRight;
end

if hrLeft <= 9 && hrLeft > 3            
    hr_ang1 = (abs((hrLeft - 9)) * (pi/6)) + pi;
elseif hrLeft <= 3 
    hr_ang1 = (3 - hrLeft) * (pi/6);
elseif hrLeft > 9
    hr_ang1 = pi - ((hrLeft - 9) * (pi/6));
end
thetaLeft = hr_ang1;     % angle in radians

scale = 15;
rho = dist * scale;     % Distance from origin to ROI
[dxRight,dyRight] = pol2cart(thetaRight, rho); 
[dxLeft,dyLeft] = pol2cart(thetaLeft, rho); % Convert tumor location as angle & dist to pixel location  
sqSize = sqSize * scale;


for j = 1:15
    xRight = xR{j} + dxRight;                       % Coordinates of center of tumor
    yRight = yR{j} - dyRight; 
    xLeft = xL{j} + dxLeft;
    yLeft = yL{j} - dyLeft;
       
    rightMinX{j} = round(xRight - sqSize/2);
    rightMinY{j} = round(yRight - sqSize/2);
    
    leftMinX{j} = round(xLeft - sqSize/2);
    leftMinY{j} = round(yLeft - sqSize/2);
    
    c1R{j} = [round(xRight - sqSize/2), round(yRight - sqSize/2)];  % Top Left Corner
    c2R{j} = [round(xRight - sqSize/2), round(yRight + sqSize/2)];  % Bottom Left Corner
    c3R{j} = [round(xRight + sqSize/2), round(yRight + sqSize/2)];  % Bottum Right Corner
    c4R{j} = [round(xRight + sqSize/2), round(yRight - sqSize/2)];  % Top Right Corner      
    
    c1L{j} = [round(xLeft - sqSize/2), round(yLeft - sqSize/2)];  % Top Left Corner
    c2L{j} = [round(xLeft - sqSize/2), round(yLeft + sqSize/2)];  % Bottom Left Corner
    c3L{j} = [round(xLeft + sqSize/2), round(yLeft + sqSize/2)];  % Bottum Right Corner
    c4L{j} = [round(xLeft + sqSize/2), round(yLeft - sqSize/2)];  % Top Right Corner      
    
    rightROI{j} = imcrop(I_mat{j}, [rightMinX{j} rightMinY{j} sqSize sqSize]);
    leftROI{j} = imcrop(I_mat{j}, [leftMinX{j} leftMinY{j} sqSize sqSize]);
    
end 

%% Find data for ROI 

for k = 1:15
    nonzeroR = rightROI{k}(find(rightROI{k}>0));
    nonzeroL = leftROI{k}(find(leftROI{k}>0));
    aveRight(k) = mean2(nonzeroR);
    aveLeft(k) = mean2(nonzeroL);
end 

for l = 1:14
    stepChangeRight(l) = aveRight(l+1) - aveRight(l);
    stepChangeLeft(l) = aveLeft(l+1) - aveLeft(l);
end


totalChangeRight = aveRight(15) - aveRight(1);
aveStepChangeRight = mean2(stepChangeRight);

totalChangeLeft = aveLeft(15) - aveLeft(1);
aveStepChangeLeft = mean2(stepChangeLeft);

%Print Results
fprintf('Clock Hour (Right): %d \nDistance: %d \nSquare Size: %d \n'...
    ,hrRight, dist, sqSize/15);
fprintf('Total Change Right: %f \nTotal Change Left: %f \nAverage Rate of Change Right: %f \nAverqage Rate of Change Left: %f\n',...
    totalChangeRight, totalChangeLeft, aveStepChangeRight, aveStepChangeLeft);

%% Verify the regions
%{
for m = 1:15
    figure, 
    I1 = I_mat{m}(find(I_mat{m}>0));
    imshow(I_mat{m}, [min(I1) max(I1)])
    hold on,
    plot(xR{m}, yR{m}, '+');
    hold on,
    plot(xL{m}, yL{m}, '+');
    hold on
    plot([c1R{m}(1 ) c2R{m}(1)],[c1R{m}(2) c2R{m}(2)],'b');                      % Create red box region on Image Display
    plot([c2R{m}(1) c3R{m}(1)],[c2R{m}(2) c3R{m}(2)],'b');
    plot([c3R{m}(1) c4R{m}(1)],[c3R{m}(2) c4R{m}(2)],'b');
    plot([c4R{m}(1) c1R{m}(1)],[c4R{m}(2) c1R{m}(2)],'b');
    hold on,
    plot([c1L{m}(1 ) c2L{m}(1)],[c1L{m}(2) c2L{m}(2)],'r');                      % Create red box region on Image Display
    plot([c2L{m}(1) c3L{m}(1)],[c2L{m}(2) c3L{m}(2)],'r');
    plot([c3L{m}(1) c4L{m}(1)],[c3L{m}(2) c4L{m}(2)],'r');
    plot([c4L{m}(1) c1L{m}(1)],[c4L{m}(2) c1L{m}(2)],'r');
end
%}