%% Loading JPEG Image

imID1 = input('Enter image name you want to open: ','s');
imID = strcat(imID1,'.jpg');

dir = uigetdir;
J = imread([dir '\' imID]);

figure, imshow(J);
title('Visible Image');

% Reshaping

B = rgb2gray(J);
figure, imshow(B);
title('Visible TIFF Image');

B = imresize(B, [512 640]); 
figure, imshow(B);
title('Resized Image');

filename = strcat(imID1, '.tif');
imwrite(B,filename);

%% Registration of JPEG 

viss = input('Enter image name you want to open: ','s'); 
vis = strcat(viss,'.tif'); 

dir1 = uigetdir; 
Ivis = imread([dir1 '\' vis]); 

reff = input('Enter image name you want to open: ','s'); 
ref = strcat(reff,'.tif'); 

dir2 = uigetdir; 
Iref = imread([dir2 '\' ref]); 

fprintf('Pick right nipple. \n');
figure, imshow(Ivis, []), title('Right Nipple Visible')
[X1,Y1] = ginput(1);

fprintf('Pick right nipple. \n');
figure, imshow(Iref, []), title('Right Nipple Reference')
[X2,Y2] = ginput(1);

val = mean(mean(Ivis));

% Cropping

hshift = round(X2 - X1);
if hshift>0
    Ivis = padarray(Ivis, [0 abs(hshift)], val, 'pre');
    Ivis = imcrop(Ivis, [0 0 640 512]);
else
    Ivis = padarray(Ivis, [0 abs(hshift)], val, 'post');
    [r c] = size(Ivis);
    Ivis = imcrop(Ivis, [(c-640) (512-r) (c) (r)]);
end

vshift = round(Y2 - Y1);
if vshift>0
    Ivis = padarray(Ivis, abs(vshift), val, 'pre');
    Ivis = imcrop(Ivis, [0 0 640 512]);
else
    Ivis = padarray(Ivis, abs(vshift), val, 'post');
    [r c] = size(Ivis);
    Ivis = imcrop(Ivis, [(640-c) (r-511) c r]); 
end

figure(1), imshow(Ivis)
figure(2), imshow(Iref, [])

% Scaling Down
% 
% fprintf('Pick right nipple. \n');
% figure, imshow(Ivis, []), title('Right Nipple Visible')
% [Xr1,Yr1] = ginput(1);
% 
% fprintf('Pick right nipple. \n');
% figure, imshow(Iref, []), title('Right Nipple Reference')
% [Xr2,Yr2] = ginput(1);
% 
% fprintf('Pick left nipple. \n');
% figure, imshow(Ivis, []), title('Left Nipple Visible')
% [Xl1,Yl1] = ginput(1);
% 
% fprintf('Pick left nipple. \n');
% figure, imshow(Iref, []), title('Left Nipple Reference')
% [Xl2,Yl2] = ginput(1);
%  

% visdist = sqrt((Xr1^2)-(Xl1^2));
% refdist = sqrt((Xr2^2)-(Xl2^2));
% 
% if visdist > refdist
%     ratio = refdist/visdist;
%     Ivis = imresize(Ivis, ratio);
% else 
%     ratio = refdist/visdist;
%     Ivis = imresize(Ivis, ratio);
%     Ivis = imcrop(Ivis, [0 0 640 512]);
% end 


% Compare

II = imshowpair(Iref, Ivis)
%% 


% Manual Segmentation 

imshow(Ivis, []);
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

blackMaskedImageR = Iref;
blackMaskedImageR(~binaryImage) = 0;

% Now crop the image.

leftColumn = min(x);
rightColumn = max(x);
topLine = min(y);
bottomLine = max(y);
width = rightColumn - leftColumn + 1;
height = bottomLine - topLine + 1;    

%newCropV = imcrop(blackMaskedImageV, [leftColumn, topLine, width, height]);
close;
newCropR = imcrop(blackMaskedImageR, [leftColumn, topLine, width, height]);
close; 

% save

newfile = strcat(reff, '_reg.tif');
imwrite(blackMaskedImageR, sprintf(newfile));  



