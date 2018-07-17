path = uigetfile('.tif');

I = imread(path);

I1 = I(find(I>0));
imshow(I, []);
hold on
hFH = imfreehand();
binaryImage = hFH.createMask();
xy = hFH.getPosition;

% Get coordinates of the boundary of the freehand drawn region.
structBoundaries = bwboundaries(binaryImage);
xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
x = xy(:, 2);   % Columns.
y = xy(:, 1);   % Rows.
% Mask the image outside the mask, and display it.

% Will keep only the part of the image that's inside the mask, zero outside mask.
blackMaskedImage = I;
blackMaskedImage(~binaryImage) = 0;

% Now crop the image.
leftColumn = min(x);
rightColumn = max(x);
topLine = min(y);
bottomLine = max(y);
width = rightColumn - leftColumn + 1;
height = bottomLine - topLine + 1;    

newCrop = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);
close;

imwrite(blackMaskedImage, sprintf('Cropped_1799_P10.tif'));  