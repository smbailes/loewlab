% I = imread('P12_CroppedContrast.tif');
% imshow(I)
% title('Original Image')
% mask = zeros(size(I));
% mask(25:end-25,25:end-25) = 1;
% figure
% imshow(mask)
% title('Initial Contour Location')
% bw = activecontour(I,mask,300);
% figure
% imshow(bw)
% title('Segmented Image')

I = imread('coins.jpg');
imshow(I)
title('Original Image')
grysc=rgb2gray(I);
figure
imshow(grysc)
title('Grayscale of Original Image')
J = imcomplement(grysc);
figure
imshow(J);
title('Complement')
mask = zeros(size(grysc));
mask(10:end-10,10:end-10) = 1;
figure
imshow(mask)
title('Initial Contour Location')
bw = activecontour(grysc,mask,10);
figure
imshow(bw)
title('Segmented Image')