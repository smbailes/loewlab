% function newcrop = RegisteredCropScript(location)

answer = questdlg('ID Patient Type:','Patient Type','Patient','Volunteer','Patient');

if (strcmp(answer, 'Patient'))
    ptID = patientselect;    % Dialog Box for patient selection
    dir = uigetdir;
    location = strcat(dir, '\', ptID, '\Cropped Image\');
else
    ptID = volunteerselect;
    dir = uigetdir;
    location = strcat(dir, '\', ptID, '\Cropped Image\');
end 

a=inputdlg('Enter first image number: '); 
strt=a{1};  
strt = str2double(strt); 
for i=1:15          
    I_mat{i} = imread([location sprintf('%04d.tif',strt)]);    % Read each image into I_mat
    strt=strt+120;            % Go to next image (for cropped), HAS TO BE CHANGED TO INCREMENT BY 1
end

%% Apply Crop to All Registered Images based on User Drawn Input
% cd ([homedir 'Registered/' ptID 'Registered/']);
firstImage = imread(location);
figure
set(gcf,'units','inches', 'Position',[4 2 10 8])
imshow(firstImage,[]);
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
blackMaskedImage = firstImage;
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

mkdir ALGCropped;
cd ALGCropped;

imwrite(newCrop,'0000.tif');

for i = 120:120:1680
    cd ..
    newImage = imread(sprintf('%04d.tif',i));
   
    blackMaskedImage = newImage;
    blackMaskedImage(~binaryImage) = 0;
    newCrop = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);

    cd ALGCropped    
    imwrite(newCrop, sprintf('%04d.tif',i));    
end

