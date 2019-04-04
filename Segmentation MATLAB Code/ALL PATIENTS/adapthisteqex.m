I=imread("0000 - P11.tif");
J=adapthisteq(I);
title('Just adapthisteq(I)');
imshow(J,[])
K=adapthisteq(I,0.02);
figure
imshow(K)
