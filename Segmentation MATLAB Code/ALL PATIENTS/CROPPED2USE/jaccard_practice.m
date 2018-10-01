 close all
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
J = imread('SnakesP16C.tif');
[x,y,z]=size(J) %necessary for snakes to run
J(:,:,1)=[]; %neccessary for snakes to run
%J(:,:,2)=[]; %neccessary for snakes to run
%J(:,:,3)=[]; %neccessary for snakes to run
%J(:,:,4)=[]; %neccessary for snakes to run
% J(:,641,:)=[]; %necessary if snakes file is bigger than 640

JJ=logical(J);
% [a,b,c]=size(JJ);
JJ(:,:,1)=[];
JJ(:,:,2)=[];
figure
%imshow(JJ)

title('Snakes') %or connectPixels
%mask = false(size(J));
%mask(25:end-25,25:end-25) = true;
%BW = activecontour(J, mask, 300);
%figure
%imshow(BW)
%title('Active Contour')
M = imread('Manual_0000_P16C.tif');
BW_groundTruth=logical(M)
figure
imshow(BW_groundTruth);
title('Truth')
similarity = jaccard(JJ, BW_groundTruth);
figure
imshowpair(JJ, BW_groundTruth)
j=similarity*100
fprintf('Jaccard Comparison for snakes and Truth. JIP= %f',j)
