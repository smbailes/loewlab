clc
clear

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
Img = imread([dir '\' ptID]); ;
[N,M] = size(Img);

filter = fspecial('average',100);%fspecial('gaussian', [3,3],100);
Img2 = imfilter(Img,filter);
BImg = im2double(Img2);
binary_thre = graythresh(BImg);
binary_img = im2bw(BImg,binary_thre);
Threshold = binary_thre*(2^16-1)+1000;

Img(Img<Threshold) = 0;
Img(Img==0) = min(Img(Img>0));
 normalized = max(Img(:))-Threshold;%min(Img(:));%Threshold;
 Img = Img-Threshold;%min(Img(:));%Threshold;
for n = 1:N
    for m = 1:M
       Img8b(n,m) = double(Img(n,m))/double(normalized);
    end
end
Img8b = uint8(Img8b.*255);
% Img8b(Img8b>0)=255; % for groundtruth
imshow(Img8b)

imwrite(Img8b,'1799_P7_8b.tif')


