%% Clear Images
close all;
clc;

%% Choose Image to input

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
I = imread([dir '\' ptID]);

%% Run Snakes

snakimg = snake(I);
imwrite(snakimg, 'SnakeImg_0000_P12.tif');