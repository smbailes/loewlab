% close all
% A = imread('0000 - P12.tif');
% J=rgb2gray(A)
% figure
% imshow(J)
% title('Original Image')
% mask = false(size(J));
% mask(25:end-25,25:end-25) = true;
% BW = activecontour(J, mask, 300, 'edge');
% figure
% imshow(BW)
% title('Active Contour')
% M = imread('Cropped_0000_P12.tif');
% BW_groundTruth=logical(M)
% similarity = jaccard(BW, BW_groundTruth);
% figure
% imshowpair(BW, BW_groundTruth)
% title(['Jaccard Index = ' num2str(similarity)])
J = imread('P4.tif');
figure
imshow(J)
title('connectPixels')
mask = false(size(J));
mask(25:end-25,25:end-25) = true;
BW = activecontour(J, mask, 300);
figure
imshow(BW)
title('Active Contour')
M = imread('Cropped_0000_P4.tif');
BW_groundTruth=logical(M)
similarity = jaccard(J, BW_groundTruth);
figure
imshowpair(J, BW_groundTruth)
