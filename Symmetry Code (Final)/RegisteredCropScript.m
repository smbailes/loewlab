%% AutoCropScript
% Applies uniform Crop to All Registered Images for a Pt

ptID = patientselect;% obtain first value in answer matrix
user = userselect;          % Dialog Box for user selection

a = strcmp(user,'Aidan');   % Compare user input to 'Aidan"
b = strcmp(user,'Lovelace');

if (a == 1 && b == 0)                 % If Aidan was selected
    location = (['/Users/AidanMurray/Data/GWBox/GRP_Loew-Doc/NadaKamona/Clinic Patients/' ptID '/']);
    homedir = '/Users/AidanMurray/Data/GWBox/GRP_Loew-Doc/NadaKamona/Clinic Patients/';
elseif (a == 0 && b == 0)
    location = (['\Users\shann\Box\GRP_Loew-Doc\NadaKamona\Clinic Patients\Manual Crop\' ptID 'Cropped\Cropped Image\']);
    homedir = '\Users\shann\Box\GRP_Loew-Doc\NadaKamona\Clinic Patients\';
elseif (a == 0 && b == 1)
    location = (['D:\BreastPatients\' ptID '\']);
    homedir = 'D:\BreastPatients\';
end

%% Check if Registered Images Exist
% cd(location);
cd([homedir 'Registered/'])
n = exist([ptID 'Registered'], 'dir');

if n ~= 7 %if the images havent been registered, call function
    registerEveryMinute(user,ptID);
end

%% Apply Crop to All Registered Images based on User Drawn Input
cd ([homedir 'Registered/' ptID 'Registered/']);
firstImage = imread('0000.tif');
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

