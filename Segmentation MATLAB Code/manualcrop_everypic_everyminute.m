%% Online Source for cropping an image manually
% Looks great!! Works perfectly
% Crops one image manually every minute. Have to crop 15 images per patient
% manually
% ***Extra parts from the original code were removed***

% Edited by : Nada Kamona
% Last Modified: 03/20/2017   11:46 AM

% Edited by: Kate Fergusson
% Last Modified: 06/21/207    11:25 AM

% Link: https://www.mathworks.com/matlabcentral/answers/225485-how-can-remove-all-regions-outside-roi

%% Part 1: Crop one image manually

% Demo to have the user freehand draw an irregular shape over a gray scale image.
% Then it creates new images:
% (1) where the drawn region is all white inside the region and untouched outside the region,
% (2) where the drawn region is all black inside the region and untouched outside the region,
% (3) where the drawn region is untouched inside the region and all black outside the region.
% It also (4) calculates the mean intensity value and standard deviation of the image within that shape,
% (5) calculates the perimeter, centroid, and center of mass (weighted centroid), and
% (6) crops the drawn region to a new, smaller separate image.

% Change the current folder to the folder of this m-file.
if(~isdeployed)
	cd(fileparts(which(mfilename)));
end
clc;	% Clear command window.
clear;	% Delete all variables.
% close all;	% Close all figure windows except those created by imtool.
% imtool close all;	% Close all figure windows created by imtool.
workspace;	% Make sure the workspace panel is showing.
fontSize = 16;

ptID = input('Enter patient ID number: ','s');
ptNum = ptID(end);

i = 0; %image counter, pics one image every minute
for n=0:14
    if i==0
        grayImage = imread(['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Patients\' ptID '\0000.tif']);
    elseif (i<1000 && i>0)
        nextImg = imread(['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Patients\' ptID '\0' num2str(i) '.tif']);  
    else 
        nextImg = imread(['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Patients\' ptID '\' num2str(i) '.tif']);
    end
    
    ImTitle = sprintf('Original Grayscale Image # %i',n);
    figure, imshow(grayImage, []), title(ImTitle, 'FontSize',fontSize);
    axis on;
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

    % Ask user to draw freehand mask.
    message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
    uiwait(msgbox(message));
    hFH = imfreehand(); % Actual line of code to do the drawing.

    % Create a binary image ("mask") from the ROI object.
    binaryImage = hFH.createMask();
    xy = hFH.getPosition;

    % Get coordinates of the boundary of the freehand drawn region.
    structBoundaries = bwboundaries(binaryImage);
    xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
    x = xy(:, 2); % Columns.
    y = xy(:, 1); % Rows.

    % Mask the image outside the mask, and display it.
    % Will keep only the part of the image that's inside the mask, zero outside mask.
    blackMaskedImage = grayImage;
    blackMaskedImage(~binaryImage) = 0;

    % Now crop the image.
    leftColumn = min(x);
    rightColumn = max(x);
    topLine = min(y);
    bottomLine = max(y);
    width = rightColumn - leftColumn + 1;
    height = bottomLine - topLine + 1;
    croppedImage = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);
    croppedOriginal = imcrop(grayImage, [leftColumn, topLine, width, height]);

    %Dispay cropped image
    figure, imshow(croppedImage), title(sprintf('Cropped Image # %i',n), ...
        'FontSize', fontSize)

        %Export/Save the image to the directory folder
    if i == 0
        imwrite(croppedImage,['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Patients\Manual Crop\' ptID 'Cropped\Cropped Image\0000.tif']);
        imwrite(croppedOriginal,['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Patients\Manual Crop\' ptID 'Cropped\Cropped Original\0000.tif']);
    elseif (i<1000 && i>0)
        imwrite(croppedImage,['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Patients\Manual Crop\' ptID 'Cropped\Cropped Image\0' num2str(i) '.tif']);
        imwrite(croppedOriginal,['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Patients\Manual Crop\' ptID 'Cropped\Cropped Original\0' num2str(i) '.tif']);
    else
        imwrite(croppedImage,['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Patients\Manual Crop\' ptID 'Cropped\Cropped Image\' num2str(i) '.tif']);
        imwrite(croppedOriginal,['\\admsrv.seas.gwu.edu\administration\Loew-Doc\NadaKamona\Clinic Patients\Manual Crop\' ptID 'Cropped\Cropped Original\' num2str(i) '.tif']);
    end
   
    i = i + 120; %increment to go to the next minute
    close all;
end

disp('Wow, wasn''t that the most tedious task in the world?')