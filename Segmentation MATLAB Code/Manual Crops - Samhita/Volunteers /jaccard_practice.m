close all
J = imread('SnakesV19C.tif');
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
%figure
%imshow(JJ)

%title('Snakes') %or connectPixels
%mask = false(size(J));
%mask(25:end-25,25:end-25) = true;
%BW = activecontour(J, mask, 300);
%figure
%imshow(BW)
%title('Active Contour')
M = imread('EManual_0000_V19C.tif'); %emilies
BW_groundTruth=logical(M)
%figure
%imshow(BW_groundTruth);
%title('Truth')
similarity = jaccard(JJ, BW_groundTruth);
%figure
imshowpair(JJ, BW_groundTruth)
j=similarity
%title('Jaccard Comparison for Snakes and Truth.')

