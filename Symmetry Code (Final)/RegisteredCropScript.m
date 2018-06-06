% function newcrop = RegisteredCropScript(location)

path = uigetdir;
location = strcat(path, '\');

strt = 120;
for i=1:14          
    I_mat{i} = imread([location sprintf('%04d.tif',strt)]);    % Read each image into I_mat
    strt=strt+120;            % Go to next image (for cropped), HAS TO BE CHANGED TO INCREMENT BY 1
end

%% Apply Crop to All Registered Images based on User Drawn Input
% cd ([homedir 'Registered/' ptID 'Registered/']);
image = I_mat{7};
I = getMatrixOutliers(image);
nonzero = I(find(I>0));
h = max(nonzero);
l = min(nonzero);

figure
set(gcf,'units','inches', 'Position',[4 2 10 8])
imshow(image,[l h]);
%     imcontrast
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
width = rightColumn - leftColumn + 1;
height = bottomLine - topLine + 1;    

newCrop = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);
close;

newLocation = strcat(location, '\', 'Cropped');
mkdir(newLocation);
cd(newLocation);

imwrite(newCrop,'0120.tif');

for i = 240:120:1680
    cd(location)
    newImage = imread(sprintf('%04d.tif',i));
   
    blackMaskedImage = newImage;
    blackMaskedImage(~binaryImage) = 0;
    newCrop = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);

    cd(newLocation)    
    imwrite(newCrop, sprintf('%04d.tif',i));    
end

