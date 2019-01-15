%% Setup: Find directory, open image
close all; clear all; clc;
location = uigetdir;
location = strcat(location, '\');

%% Read in Images
strt = 0;
for i=1:15          
    I_mat{i} = imread([location sprintf('%04d.tif',strt)]);    % Read each image into I_mat
    strt=strt+120;            % Go to next image (for cropped), HAS TO BE CHANGED TO INCREMENT BY 1
end

%% Show images to determine whether to use affine or demons registration
image1 = I_mat{1};
image2 = I_mat{15};
I = getMatrixOutliers(image1);
nonzero = I(find(I>0));
h = max(nonzero);
l = min(nonzero);


figure, subplot(2,1,1),
imshow(image1, [l h]);
subplot(2,1,2),
imshow(image2, [l h]);

pause
close all

%% Register Images
%Use demons for patients that move more or with larger breasts

in = input('Affine or demons? Enter a/d: ','s');
fixed = image1;


if strcmp(in, 'a')%Affine registration

    
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumStepLength = 0.04;
    optimizer.MaximumIterations = 100; %SETTING FOR NEW MATLAB

    registeredImage = imregister(image2,fixed,'affine',optimizer,metric);
    fprintf('Finished Affine Registration\n');
    
elseif strcmp(in, 'd') %demons registration   


    [~,registeredImage] = imregdemons(image2,fixed,[500 400 200],'AccumulatedFieldSmoothing',3);

    %Bring registered image back to CPU 
    fprintf('Finished Demons Registration\n'); 
end

%% Show registered images 
figure, subplot(2,1,1)
imshow(fixed, [l h]);
subplot(2,1,2)
imshow(registeredImage, [l h]);
pause
close all
%% Crop images
figure,
imshow(fixed, [l h]);

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
BW1 = fixed;
BW1(~binaryImage) = 0;

BW2 = registeredImage;
BW2(~binaryImage) = 0;

% Now crop the image.
leftColumn = min(x);
rightColumn = max(x);
topLine = min(y);
bottomLine = max(y);
width = rightColumn - leftColumn + 5;
height = bottomLine - topLine + 5;    

cropped1 = imcrop(BW1, [leftColumn, topLine, width, height]);
cropped2 = imcrop(BW2, [leftColumn, topLine, width, height]);

close;


%% Show cropped images 
figure, subplot(2,1,1)
imshow(cropped1, [l h]);
subplot(2,1,2)
imshow(cropped2, [l h]);
pause
close all

%% Write
newLocation = strcat(location, 'Cropped 1 and 15\');
mkdir(newLocation);
cd(newLocation);
imwrite(cropped1, '0000.tif');
imwrite(cropped2, '1680.tif');

subtract;