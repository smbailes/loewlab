clear all;
close all;
clc;
%% Patient Selection
    [location, ptID] = pathfinder; 
    newLocation = strcat(location, '\', 'Cropped');
    mkdir(newLocation);
    
%% Read in Images
    cd(location);
    left0 = imread([ptID '_Left_min0.tif']);
    right0 = imread([ptID '_Right_min0.tif']);
    left15 = imread([ptID '_Left_min15.tif']);
    right15 = imread([ptID '_Right_min15.tif']);
    
    I_mat{1} = left0; I_mat{2} = right0; 
    I_mat{3} = left15; I_mat{4} = right15;

%% Crop Images
for i = 1:4
    image = I_mat{i};
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
    % hFH = imrect();
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

    newCrop{i} = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);
    close;
end 

%% Write new files
cd(newLocation);
imwrite(newCrop{1}, sprintf('%s_Left_min0.tif',ptID));
imwrite(newCrop{2}, sprintf('%s_Right_min0.tif',ptID));
imwrite(newCrop{3}, sprintf('%s_Left_min15.tif',ptID));
imwrite(newCrop{4}, sprintf('%s_Right_min15.tif',ptID));

%% Show Images
for j = 1:4
    figure
    imshow(newCrop{j}, [l h]);
end
