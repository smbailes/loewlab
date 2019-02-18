ptIDref = input('Enter image name you want to open: ','s'); 
ptIDref = strcat(ptIDref,'.tif'); 

dir1 = uigetdir; 
Iref = imread([dir1 '\' ptIDref]); 

%% 


ptIDvis = input('Enter image name you want to open: ','s'); 
ptIDvis = strcat(ptIDvis,'.tif'); 

dir2 = uigetdir; 
Ivis = imread([dir2 '\' ptIDvis]); 

figure(1), imshow(Iref,[]) 
title('Original Image')

figure(2), imshow(Ivis,[])
title('Original Image')