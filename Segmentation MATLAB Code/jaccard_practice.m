close all
A = imread('0000 - P1.tif');
J = im2uint8(A)
figure
imshow(J)
title('Original Image')
mask = false(size(J));
mask(25:end-25,25:end-25) = true;
BW = activecontour(J, mask, 300);
figure
imshow(BW)
title('Active Contour')
M = imread('Cropped_0000_P1.tif');
BW_groundTruth=logical(M)
similarity = jaccard(BW, BW_groundTruth);
figure
imshowpair(BW, BW_groundTruth)
title(['Jaccard Index = ' num2str(similarity)])