% %% Part 1 Initial Loading of Images
% 
% ptID = input('Enter image name you want to open: ','s'); %Request patient image name
% ptID = strcat(ptID,'.tif'); 
% %I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
% dir = uigetdir; 
% I = imread([dir '\' ptID]); 
% 
% figure, imshow(I,[]) %to help decide if it should be cropped or not
% title('Original Image')

%% Loading JPEG Image

imID1 = input('Enter image name you want to open: ','s');
imID = strcat(imID1,'.jpg');

dir = uigetdir;
J = imread([dir '/' imID]);

%figure, imshow(J);
%title('Visible Image');


%% Reshaping

B = rgb2gray(J);
figure, imshow(B);
title('Visible TIFF Image');

B = imresize(B, [512 640]); 
figure, imshow(B);
title('Resized Image');

filename = strcat(imID1, '.tif');
imwrite(B,filename);
