path1 = uigetfile('.tif');
I = imread(path1);

I = imshow(I, []);
hold on
[J1, rect] = imcrop(I);
%imwrite(J1, sprintf('Original_P7.tif')); 

%[r c] = size(J1);
%numrows = 480 - r;
%J1 = padarray(J1, numrows, 0, 'pre');


imwrite(J1, sprintf('0013_0000.tif')); 


close all
I = [];

%% 

path2 = uigetfile('.tif');
I = imread(path2);
J2 = imcrop(I, rect);
%J2 = padarray(J2, numrows, 0, 'pre');

imwrite(J2, sprintf('0013_0900.tif')); 

%% 


imwrite(J2, sprintf('Manual3_P1C.tif')); 

path3 = uigetfile('.tif');
I = imread(path3);
J3 = imcrop(I, rect);
%J3 = padarray(J3, numrows, 0, 'pre');


imwrite(J3, sprintf('0013_1799.tif')); 
%% 

figure(1)
imshow(J1, []);
figure(2)
imshow(J2, []);