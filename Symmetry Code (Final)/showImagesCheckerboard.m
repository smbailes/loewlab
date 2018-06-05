%% Clean Up
clear all;
close all;
clc; 

%% Choose reference image

newLocation = uigetdir;
% figure
%Find max and min for contrast 
path = [newLocation '\' sprintf('%04d.tif',840)];
ref = imread(path);
%% Image Contrast
% Not used with imshowpair

I = getMatrixOutliers(ref);
I_nonzero = I(find(I>0));
h = max(I_nonzero);
l = min(I_nonzero);
%% Path to Images normally registered and registered using demons registration

registeredLocation = strcat(newLocation, '\Registered\');
demonRegLocation = strcat(newLocation, '\Demons Registered\');

%% Display the Images
% imshowpair used to checkerboard the reference image and registered image
% top image is normally reigstered, bottom demons registered

for k = 0000:120:1680
    RegPath = [registeredLocation '\' sprintf('%04d.tif',k)];
    registered = imread(RegPath);
    
    DRPath = [demonRegLocation '\' sprintf('%04d.tif',k)];
    dRegistered = imread(DRPath);
    
    figure %new figure for each minute
    
    %normal registration on top
    subplot(2,1,1)
    imshowpair(registered,ref,'checkerboard','Scaling','joint'); 
    %demons registered on bottom 
    subplot(2,1,2)
    imshowpair(dRegistered,ref,'checkerboard','Scaling','joint');
end

