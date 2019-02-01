%% Conversion from 16 Bit to 8 Bit Images
close all;

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif');
dir = uigetdir; 
I = imread([dir '\' ptID]);

figure;
imshow(I, [])
pause;

I(I < 8000) = 0;
I(I == 0) = min(I(I > 0));
J = (I - 8000);
J = double(J/.2000);
J = J*.255;
imwrite(J, sprintf('ConvertBit_1799_V1.tif'));
%% Manual Crop of 8 Bit Images

K = J(find(J > 0));
figure, imshow(J, [])
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

imwrite(blackMaskedImage, sprintf('ManualCrop8Bit_0000_P1.tif')); 