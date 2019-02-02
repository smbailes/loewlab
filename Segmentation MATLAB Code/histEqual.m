
ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
dir = uigetdir; 
I = imread([dir '/' ptID]); 
%%I = im2int8(I1); 

J = adapthisteq(I,'clipLimit',0.02,'Distribution','rayleigh');
%higher clip limit = more contrast 

%J = adapthisteq(I);
%CLAHE --> contrast limiting, prevent overamplification of noise

imwrite(J, 'image.png'); 
figure(1)
imshowpair(I,J,'montage')
axis off

figure(2)
imhist(I,64)

figure(3)
imhist(J,64)

