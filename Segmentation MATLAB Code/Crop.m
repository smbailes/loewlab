% ptID = input('Enter original image name you want to open: ','s'); %Request patient image name
% ptID = strcat(ptID,'.tif'); 
% %I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
% dir = uigetdir; 
% I = imread([dir '/' ptID]);
% 
% figure, imshow(I,[]);
% [x,y] = ginput(1);

mID = input('Enter M image: ','s'); %Request patient image name
mID = strcat(mID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
G = imread([dir '/' mID]);
G;

cpID = input('Enter CP image: ','s'); %Request patient image name
cpID = strcat(cpID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
J = imread([dir '/' cpID]);
J = J(1:640, y:512);

cpvID = input('Enter CPV image: ','s'); %Request patient image name
cpvID = strcat(cpvID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
K = imread([dir '/' cpvID]);
K = K(1:640, y:512);

sID = input('Enter S image: ','s'); %Request patient image name
sID = strcat(sID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
L = imread([dir '/' sID]);
L = L(1:640, y:512);

svID = input('Enter SV image: ','s'); %Request patient image name
svID = strcat(svID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
N = imread([dir '/' svID]);
N = N(1:640, y:512);