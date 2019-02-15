%% Part 1 Intial Loading of Images

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
I = imread([dir '\' ptID]); 

figure, imshow(I,[]) %to help decide if it should be cropped or not
title('Original Image')

%% Loading JPEG Image

imID = input('Enter image name you want to open: ','s');
imID = strcat(imID,'.jpg');

dir = uigetdir;
J = imread([dir '\' imID]);

figure, imshow(J);
title('Visible Image');

%% Reshaping

B = reshape(J,480,640);

