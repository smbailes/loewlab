close all;

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif');
dir = uigetdir; 
I = imread([dir '\' ptID]);

figure;
imshow(I, [])
pause;

I(I < 8000) = 0;
I(I == 0) = min(I(I > 0));
J = (I - 8000);
J = double(J/.2000);
J = J*.255;

figure;
imshow(J, [])
imwrite(J, sprintf('ConvertBit_0000_P1.tif'));