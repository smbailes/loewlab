 close all
M1 = imread('EManual_0000_P12C.tif'); %emilies crop
BW1_groundTruth=logical(M1)

M2 = imread('0000 - P12C.tif'); %samhitas crop
BW2_groundTruth=logical(M2)
%figure
%imshow(BW2_groundTruth);
%title('Truth')
similarity = jaccard(BW1_groundTruth, BW2_groundTruth);
%figure
imshowpair(BW1_groundTruth, BW2_groundTruth)
j=similarity
%title('Jaccard Comparison for Snakes and Truth')
