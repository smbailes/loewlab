ptID = input('Enter image name you want to open: ','s'); % ENTER PATIENT IMAGE NUMBER
%ex. 0000, 0120, 0240, 0360, 0480...1800 or 1799
ptID = strcat(ptID,'.tif'); % ptID labelled as TIF file (ex. 0000.tif)
dir = uigetdir; % select patient number's folder from folder with 8 bit images
I = imread([dir '\' ptID]); % retrieves image from selected folder 
I1 = I(find(I>0));

imshow(I, []); % displays image 
hold on
hFH = imfreehand(); % DRAW OUTLINE TO MANUALLY SEGMENT
% WHEN TRACING:
    % 1. start from the armpit on the left side of the screen (the patient's
    % right) (hold down the mouse and drag)
    % 2. trace the outline of the underside of the breast 
    % 3. end at the armpit on the right side of the screen (the patient's
    % left) (release the mouse it saves automatically)
    
binaryImage = hFH.createMask(); % creates binary mask based on outline
xy = hFH.getPosition; % gets boundaries of outline

% Get coordinates of the boundary of the freehand drawn region.
structBoundaries = bwboundaries(binaryImage); % Get coordinates of the boundary of outline
xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
x = xy(:, 2);   % Columns.
y = xy(:, 1);   % Rows.
% Mask the image outside the mask, and display it.

% Will keep only the part of the image that's inside the mask, zero outside mask.
blackMaskedImage = I; %makes blackMaskedImage thats equal to the original image
blackMaskedImage(~binaryImage) = 0; %this sets what was not in the traced region to 0 

% Now crop the image.
leftColumn = min(x); % this crops the image so that the left column in the minimum x value (for the traced region)
rightColumn = max(x); % this crops the image so that the right column in the maximum x value (for the traced region)
topLine = min(y); % this crops the image so that the top row is the minimum y value (for the traced region)
bottomLine = max(y); % this crops the image so that the top row is the maximum y value (for the traced region)
width = rightColumn - leftColumn + 1; % sets the width for the traced region 
height = bottomLine - topLine + 1;  % sets the height for the traced region   

newCrop = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]); %this crops teh blackmasked image to this size 
close; 

imwrite(blackMaskedImage, sprintf('IRVT029_0015M.tif')); % CHANGE THIS BEFORE RUNNING CODE
% format: 'PatientNumber_ImageNumberM.tif' 
% ex. for Patient Number 1 and Image Number 0000: 'P1_0000M.tif'