path1 = uigetfile('.tif');
I = imread(path1);

I = imshow(I, []);
hold on
[J1, rect] = imcrop(I);
imwrite(J1, sprintf('Original3_P1C.tif')); 

close all
I = [];

%% 

path2 = uigetfile('.tif');
I = imread(path2);
J2 = imcrop(I, rect);

imwrite(J2, sprintf('Manual3_P1C.tif')); 

%% 

figure(1)
imshow(J1, []);
figure(2)
imshow(J2, []);