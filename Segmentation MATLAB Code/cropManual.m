path1 = uigetfile('.tif');
I = imread(path1);

I = imshow(I, []);
hold on
[J1, rect] = imcrop(I);
imwrite(J1, sprintf('ex1.tif')); 

close all
I = [];

path2 = uigetfile('.tif');
I = imread(path2);
I = imshow(I, []);
J2 = imcrop(I, rect);

imwrite(J2, sprintf('ex2.tif')); 



figure(1)
imshow(J1, []);
figure(2)
imshow(J2, []);