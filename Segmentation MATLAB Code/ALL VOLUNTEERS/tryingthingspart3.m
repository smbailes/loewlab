snakeID = input('Enter image name you want to open: ','s'); %Request patient image name
snakeID = strcat(snakeID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Snakes Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
I = imread([dir '/' snakeID]); 

figure, imshow(I,[]) %to help decide if it should be cropped or not
[r,c]=size(I)

[x,y] = getpts
%locminy=y(find(min(y)))
miny=min(y)
% c(:,0:locminy)=0
%  c=0
% while c<locminy
%     y=0;
%     c+1;
% end
imshow(I(miny:480,:), [])
