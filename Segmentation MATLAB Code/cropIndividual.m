ptID = input('Enter image name you want to open: ','s'); %Request patient image name 
    % ^ this will be 0000, 0120, 0240, 0360, 0480, 0600, 0720, 0840, 0960,
    % 1080, 1200, 1320, 1440, 1560, 1680, 1800 or 1799
ptID = strcat(ptID,'.tif'); %makes ptID the number you enter in "dot" tif (ex: 0000.tif)
dir = uigetdir; %you will go into the folder with 8 bit images, then the patient number's folder and click open
I = imread([dir '/' ptID]); %reads in the image from the folder you selected
I1 = I(find(I>0));
imshow(I, []); %shows the image 
hold on
hFH = imfreehand(); %this is when you draw the outline to crop
%WHEN TRACING:
    %1. start from the armpit on the left side of the screen (the patient's
    %right) (hold down the mouse and drag)
    %2. trace the outline of the underside of the breast 
    %3. end at the armpit on the right side of the screen (the patient's
    %left) (release the mouse it saves automatically
binaryImage = hFH.createMask(); %creates a binary mask based on what you traced
xy = hFH.getPosition; %gets the boundaries of what you traced

% Get coordinates of the boundary of the freehand drawn region.
structBoundaries = bwboundaries(binaryImage); % Get coordinates of the boundary of the tracing
xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
x = xy(:, 2);   % Columns.
y = xy(:, 1);   % Rows.
% Mask the image outside the mask, and display it.

% Will keep only the part of the image that's inside the mask, zero outside mask.
blackMaskedImage = I; %makes blackMaskedImage thats equal to the original image
blackMaskedImage(~binaryImage) = 0; %this sets what was not in the traced region to 0 

% Now crop the image.
leftColumn = min(x); %this crops the image so that the left column in the minimum x value (for the traced region)
rightColumn = max(x); %this crops the image so that the right column in the maximum x value (for the traced region)
topLine = min(y); %this crops the image so that the top row is the minimum y value (for the traced region)
bottomLine = max(y); %this crops the image so that the top row is the maximum y value (for the traced region)
width = rightColumn - leftColumn + 1; %sets the width for the traced region 
height = bottomLine - topLine + 1;  %sets the height for the traced region   

newCrop = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]); %this crops teh blackmasked image to this size 
close; 



imwrite(blackMaskedImage, sprintf('hi_Manual3_P11.tif'));  %CHANGE THIS
% ^ change to PatientNumber_ImageNumberM.tif (ex: for Patient 1 0000:
% P1_0000M