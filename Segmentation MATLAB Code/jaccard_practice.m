close all
A = imread('hands1.jpg');
I = rgb2gray(A);
figure
imshow(I)
title('Original Image')
mask = false(size(I));
mask(25:end-25,25:end-25) = true;
BW = activecontour(I, mask, 300);
figure
title('Active Contour')
imshow(BW)
BW_groundTruth = imread('hands1-mask.png');
similarity = jaccard(BW, BW_groundTruth);
figure
imshowpair(BW, BW_groundTruth)
title(['Jaccard Index = ' num2str(similarity)])