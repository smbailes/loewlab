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
    
%% Test fibermetric (compare fibermetric to original image)

n = 8;
image = I_mat{n};
I = image(find(image>0));
for i = 1:20
    V = fibermetric(I_mat{n},i,'ObjectPolarity','bright','StructureSensitivity',12);
    figure, 
    subplot(2,1,1);
    imshow(I_mat{n},[min(I) max(I)]);
    subplot(2,1,2);
    imshow(V)
    xlabel(sprintf('Thickness = %f',i))
end


thickness = input('Best Thickness?');
close all;

for j = 1:30
    V = fibermetric(I_mat{n},thickness,'ObjectPolarity','bright','StructureSensitivity',j);
    figure, 
    subplot(2,1,1);
    imshow(I_mat{n},[min(I) max(I)]);
    subplot(2,1,2);
    imshow(V)
    xlabel(sprintf('Structure Sensitivity = %f',j))
end
% %% Show all fibermetrics on one subplot
% 
% n = 8;
% image = I_mat{n};
% I = image(find(image>0));
% figure
% for i = 1:20
%     V = fibermetric(I_mat{n},i,'ObjectPolarity','bright','StructureSensitivity',12);
%     subplot(5,4,i);
%     imshow(V)
% end
ss = input('Best Structure Sensitivity?');

close all
figure, 
subplot(4,1,1),
imshow(image, [min(I) max(I)]);
title('Original')
V = fibermetric(image,thickness,'ObjectPolarity','bright','StructureSensitivity',ss);
subplot(4,1,2)
imshow(V)
title('Without CLAHE')


