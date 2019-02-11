close all;

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif');
dir = uigetdir; 
I = imread([dir '\' ptID]);

figure;
imshow(I)

I(I < 8000) = 0;
I(I == 0) = min(I(I > 0));
img2 = (I - 8000);
img2 = float(img2/.2000);
img2 = int(img2*255);
figure;
imshow(img2)
Img2 = imread('D:\Breast Region by DL\1234\ORG_0.png');
Img2 = rgb2gray(Img2);

figure;
imshow(J)