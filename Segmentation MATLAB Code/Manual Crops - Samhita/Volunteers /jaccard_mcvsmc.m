 close all
M1 = imread('EManual_0000_V20C.tif'); %emilie crop
BW1_groundTruth=logical(M1)

M2 = imread('Manual_0000_V20CS.tif'); %samhita crop
BW2_groundTruth=logical(M2)
%figure
%imshow(BW2_groundTruth);
%title('Truth')
similarity = jaccard(BW1_groundTruth, BW2_groundTruth);
%figure
%imshowpair(BW1_groundTruth, BW2_groundTruth)
j=similarity
%title('Jaccard Comparison for Snakes and Truth')
