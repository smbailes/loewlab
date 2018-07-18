%Compare two sides

%% Inputs 

[location, ptID] = pathfinder; 

a=0;
for i=1:15          
    I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
    a=a+120;            % Go to next image (for cropped)
end

prompt = {'Clock Hour', 'Distance', 'Square size'};  
dlg_title = 'Select a Region';                                         % box title
defaultans = {'12','1','2'};          % default inputs
options.Resize = 'on';                                                  % allows for resizing of box
answer = inputdlg(prompt, dlg_title, [1 50], defaultans, options);      % creates box

hr = str2double(answer{1});
dist = str2double(answer{2});
sqSize = str2double(answer{3});



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
questdlg('Switch sides','Switch sides','Ok','Sure','Ok')

figure('Name','Select nipple (Left)'), 
for i =1:15
    I1 = I_mat{i}(find(I_mat{i}>0));
    imshow(I_mat{i}, [min(I1) max(I1)]) % gets coordinates of nipple
    hold on,
    xlabel('-->')
    [xL{i},yL{i}] = ginput(1);
end
close

%% Create boxes
if hr <= 9 && hr > 3            
    hr_ang = (abs((hr - 9)) * (pi/6)) + pi;
elseif hr <= 3 
    hr_ang = (3 - hr) * (pi/6);
elseif hr > 9
    hr_ang = pi - ((hr - 9) * (pi/6));
end
theta = hr_ang;     % angle in radians

scale = 15;
rho = dist * scale;     % Distance from origin to ROI
[dx,dy] = pol2cart(theta, rho); % Convert tumor location as angle & dist to pixel location  

for j = 1:15
    xRight = xR{j} + dx;                       % Coordinates of center of tumor
    yRight = yR{j} - dy; 
    xLeft = xL{j} + dx;
    yLeft = yL{j} - dy;
    
    sqSize = sqSize * scale;                    % Convert X and Y dimensions of tumor from cm to pixels
   
    rightMinX = round(xRight - sqSize/2);
    rightMinY = round(yRight - sqSize/2);
    
    leftMinX = round(xLeft - sqSize/2);
    leftMinY = round(yLeft - sqSize/2);
    
    rightROI{j} = imcrop(I_mat{j}, [rightMinX rightMinY sqSize sqSize]);
    leftROI{j} = imcrop(I_mat{j}, [leftMinX leftMinY sqSize sqSize]);
    
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

totalChangeRight = aveRight(15) - aveRight(1)
aveStepChangeRight = mean2(stepChangeRight)

totalChangeLeft = aveLeft(15) - aveLeft(1)
aveStepChangeLeft = mean2(stepChangeLeft)

    
    

