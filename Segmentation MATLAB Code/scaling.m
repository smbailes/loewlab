ptIDref = input('Enter image name you want to open: ','s'); 
ptIDref = strcat(ptIDref,'.tif'); 

dir1 = uigetdir; 
Iref = imread([dir1 '\' ptIDref]); 

%% 

ptID = input('Enter patient: ','s'); 

for i = 1000:1000:4000
    
    ptIDvis = strcat(ptID, '-', num2str(i), '.jpg'); 

    dir2 = uigetdir; 
    I = imread([dir2 '\' ptIDvis]); 

    Ivis = rgb2gray(I);
    filename = strcat('IRVT025-', num2str(i), '.tif');
    imwrite(Ivis, filename);
end

% figure(1), imshow(Iref,[]) 
% title('Original Image')
% 
% figure(2), imshow(Ivis,[])
% title('Original Image')