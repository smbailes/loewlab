clear all, close all
% %% Patient Selection
%     [location, ptID] = pathfinder; 

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread([sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end
 image = I_mat{8};    
 I = image(find(image>0)); 
 figure, 
 imshow(image, [min(I) max(I)]); 
 
%% Test fibermetric (compare fibermetric to original image)

n = 8;
image = I_mat{n};
I = image(find(image>0));
% Loop through thickness values from 1-20
% Look at results and pick the best value 
for i = 1:20
    V = fibermetric(I_mat{n},i,'ObjectPolarity','bright','StructureSensitivity',12);
    figure, 
    subplot(2,1,1);
    imshow(I_mat{n},[min(I) max(I)]);
    subplot(2,1,2);
    imshow(V)
    xlabel(sprintf('Thickness = %f',i))
end

% Which value of thickness is best? Comes from looking at images and is
% input into the next loop that checks for the best structure sensitivity 
thickness = input('Best Thickness?');
close all;

% Loop through structure sensitivity values from 1-30 
% Look at results and select the best value 
for j = 1:30
    V = fibermetric(I_mat{n},thickness,'ObjectPolarity','bright','StructureSensitivity',j);
    figure, 
    subplot(2,1,1);
    imshow(I_mat{n},[min(I) max(I)]);
    subplot(2,1,2);
    imshow(V)
    xlabel(sprintf('Structure Sensitivity = %f',j))
end

ss = input('Best Structure Sensitivity?');
close all

%% Show results with best thickness and structure sensitivities  
figure, 
subplot(4,1,1),
imshow(image, [min(I) max(I)]);
title('Original')
V = fibermetric(image,thickness,'ObjectPolarity','bright','StructureSensitivity',ss);
subplot(4,1,2)
imshow(V)
title('Without CLAHE')

%% Threshold image to find vessels 
V1 = fibermetric(image, 12, 'ObjectPolarity', 'bright', 'StructureSensitivity', 12);
figure, imshow(V1)

thresh = V1 >= 0.2;
figure, imshow(thresh)




