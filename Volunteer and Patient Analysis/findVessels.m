clear all, close all
%% Patient Selection
    [location, ptID] = pathfinder; 

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end

image = I_mat{8};
I = image(find(image>0));
figure, imshow(image,[min(I) max(I)])
%% Try CLAHE
CL = 0.02;
J = adapthisteq(I_mat{8},'ClipLimit',CL,'NBins',double((max(max(image)))));
figure
imshow(J,[min(min(J)) max(max(J))])
title(sprintf('CLAHE Results %f',CL))

%% Find vessels
V1 = fibermetric(J,'ObjectPolarity','bright');
figure
imshow(V1)

V2 = V1>=0.2;
imshow(V2)
title('Vesselness with CLAHE')

