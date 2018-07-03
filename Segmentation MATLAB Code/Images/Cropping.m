ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
I = imread([dir '\' ptID]); 

figure, imshow(I,[])

%stop ellipses from going into the bottom 30 percent of an image
%set upper and lower bounds for the cutoff of ellipses
lower = 85; upper = 15;

%make and image and make the top and bottom fifteen percent zero values

%find height and length of the image in pixels
imlen = length(I);
totpix = numel(I);
imheight = totpix/imlen;

J = I;
J(1:upper*imheight,:) = zeros;
J(lower*imheight:imheight,:) = zeros;

%run ellipseDetection on J and then put ellipses into I
%maybe use bins to stop it picking up dark points outside the body for the
%ellipses