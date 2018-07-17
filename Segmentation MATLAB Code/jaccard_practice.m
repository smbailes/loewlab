close all
A = imread('0000 - P1.tif');
J = im2uint8(A)
%A=imread('hands1.jpg')
%A=[x,y,z]
%z=NULL
%a=[x,y]
%I = rgb2gray(J);
figure
imshow(J)
title('Original Image')
mask = false(size(J));
mask(25:end-25,25:end-25) = true;
BW = activecontour(J, mask, 300);
figure
title('Active Contour')
imshow(BW)
%BW_groundTruth = imread('0000 - P1.tif');
%BW_groundTruth = imread('hands1-mask.png');
M = logical(J)
BW_groundTruth = M
similarity = jaccard(BW, BW_groundTruth);
figure
imshowpair(BW, BW_groundTruth)
title(['Jaccard Index = ' num2str(similarity)])