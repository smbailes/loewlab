close all;
clear all;
clc; 

[location, ptID,answer] = pathfinder; 
strt = 1;
for i=1:2:30          
    I_mat{strt} = imread([location sprintf('%s-%04d.tif',ptID,i)]);    % Read each image into I_mat
    strt = strt+1;
end
newLocation = strcat(location, 'Registered\');
mkdir(newLocation)

%% Show images to determine whether to use affine or demons registration
image = I_mat{8};
I = getMatrixOutliers(image);
nonzero = I(find(I>0));
h = max(nonzero);
l = min(nonzero);

for i = 1:15
    figure
    imshow(I_mat{i}, [l h]);
end
pause
close all
%% Register Images
%Use demons for patients that move more or with larger breasts

in = input('Affine or demons? Enter a/d: ','s');
fixed = I_mat{8};
cd(newLocation);

if strcmp(in, 'a')%Affine registration
    for i = 1:15
        [optimizer, metric] = imregconfig('monomodal');
        optimizer.MaximumStepLength = 0.04;
        optimizer.MaximumIterations = 100;
        
        registeredImage{i} = imregister(I_mat{i},fixed,'affine',optimizer,metric);
        imwrite(registeredImage{i},sprintf('%04d.tif',i));
    end 
    fprintf('Finished Affine Registration \n');
%     for i = 0:120:1680
%     newFile = [location '/' sprintf('%04d.tif',i)];
%     newImage = imread(newFile);
%     
%     [optimizer, metric] = imregconfig('monomodal');
%     optimizer.MaximumStepLength = 0.04;
%     optimizer.MaximumIterations = 100; %SETTING FOR NEW MATLAB
% 
%     registeredImage = imregister(newImage,fixed,'affine',optimizer,metric);
%     imwrite(registeredImage,sprintf('%04d.tif',i))
%     end
%     fprintf('Finished Affine Registration\n');
elseif strcmp(in, 'd') %demons registration
    
    fixedGPU = gpuArray(fixed); %Create gpuArray
    for i = 0:120:1680
        imgpath = [location '\' sprintf('%04d.tif',i)];
        moving = imread(imgpath);
        movingGPU = gpuArray(moving); %Create gpuArray

        [~,movingReg] = imregdemons(movingGPU,fixedGPU,[500 400 200],'AccumulatedFieldSmoothing',3);
        
        %Bring registered image back to CPU
        registeredImage = gather(movingReg);
        
        imwrite(registeredImage,sprintf('%04d.tif',i))
    end 
    fprintf('Finished Demons Registration\n'); 
end

%% Show newly registered images

for p = 1:15
    figure
    imshow(registeredImage{p}, [l h]);
end
pause 
close all;

%% Crop Images
figure,
imshow(registeredImage{7},[l h]);

message = sprintf('Left click and hold to begin outlining the breast region.\nSimply lift the mouse button to finish');
uiwait(msgbox(message));
hFH = imfreehand();
% Create a binary image ("mask") from the ROI object.
binaryImage = hFH.createMask();
xy = hFH.getPosition;

% Get coordinates of the boundary of the freehand drawn region.
structBoundaries = bwboundaries(binaryImage);
xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
x = xy(:, 2);   % Columns.
y = xy(:, 1);   % Rows.
% Mask the image outside the mask, and display it.

% Will keep only the part of the image that's inside the mask, zero outside mask.
blackMaskedImage = image;
blackMaskedImage(~binaryImage) = 0;

% Now crop the image.
leftColumn = min(x);
rightColumn = max(x);
topLine = min(y);
bottomLine = max(y);
width = rightColumn - leftColumn + 5;
height = bottomLine - topLine + 5;    

% newCrop = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);
% close;

newLocation = strcat(location, 'Cropped\');
mkdir(newLocation);
cd(newLocation);
for i = 1:15
    
    cd(location)
    newImage = registeredImage{i};
   
    blackMaskedImage = newImage;
    blackMaskedImage(~binaryImage) = 0;
    newCrop{i} = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);

    cd(newLocation)    
    imwrite(newCrop{i}, sprintf('%04d.tif',i));    
end