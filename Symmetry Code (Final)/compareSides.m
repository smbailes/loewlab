%Compare two sides
%after running it once for a patient comment out Location and Select
%nipples sections

clearvars -except location ptId I_mat xR yR xL yL,
close all, 

%% Location
[location, ptID] = pathfinder; 

a=0;
for i=1:15          
    I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
    a=a+120;            % Go to next image (for cropped)
end

%% Select Nipples

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

%% Parameters

prompt = {'Clock Hour', 'Distance', 'Square size'};  
dlg_title = 'Select a Region';                                         % box title
defaultans = {'12','1','2'};          % default inputs
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
    
    rightROI{j} = imcrop(I_mat{j}, [rightMinX{j} rightMinY{j} sqSize sqSize]);
    leftROI{j} = imcrop(I_mat{j}, [leftMinX{j} leftMinY{j} sqSize sqSize]);
    
end 

%% 
for m = 1:15
    figure, 
    I1 = I_mat{m}(find(I_mat{m}>0));
    imshow(I_mat{m}, [min(I1) max(I1)])
    hold on,
    plot(rightMinX{m}, rightMinY{m}, '*');
    hold on,
    plot(xR{m}, yR{m}, '+');
    hold on,
    plot(leftMinX{m}, leftMinY{m}, '*');
    hold on,
    plot(xL{m}, yL{m}, '+');
end

%% Find data for ROI 

for k = 1:15
    aveRight(k) = mean2(rightROI{k});
    aveLeft(k) = mean2(leftROI{k});
end 

for l = 1:14
    stepChangeRight(l) = aveRight(l+1) - aveRight(l);
    stepChangeLeft(l) = aveLeft(l+1) - aveLeft(l);
end

fprintf('Clock Hour: %d \nDistance: %d \nSquare Size: %d \n',hrRight, dist, sqSize);

totalChangeRight = aveRight(15) - aveRight(1)
aveStepChangeRight = mean2(stepChangeRight)

totalChangeLeft = aveLeft(15) - aveLeft(1)
aveStepChangeLeft = mean2(stepChangeLeft)
    
