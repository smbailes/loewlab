% I = imread('P12_CroppedContrast.tif');
% imshow(I)
% title('Original Image')
% mask = zeros(size(I));
% mask(25:end-25,25:end-25) = 1;
% figure
% imshow(mask)
% title('Initial Contour Location')
% bw = activecontour(I,mask,300);
% figure
% imshow(bw)
% title('Segmented Image')
close all
%I = imread('coins.jpg');
I=imread('P1_CroppedContrast.tif');
imshow(I)
title('Original Image')
grysc=rgb2gray(I);
%figure
%imshow(grysc)
%title('Grayscale of Original Image')
%J = imcomplement(grysc);
%figure
%imshow(J);
%title('Complement')
mask = false(size(grysc));
mask(25:end-25,25:end-25) = true;
figure
imshow(mask)
title('Initial Contour Location')
bw = activecontour(grysc,mask,300);
figure
imshow(bw)
title('Segmented Image')