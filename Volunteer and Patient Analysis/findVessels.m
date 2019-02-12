clear all, 
close all

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread(sprintf('%04d.tif',a));    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end

image = I_mat{8};
I = image(find(image>0));
%% Try CLAHE
CL = 0.025;
J = adapthisteq(I_mat{8},'ClipLimit',CL,'NBins',double((max(max(image)))));

%% Find vessels
V1 = fibermetric(J,'ObjectPolarity','bright');
V2 = V1>=0.2;

%% Compare to 8-bit 
I8 = imread('8bit.png');
img8 = I8(find(I8>0));

J8 = adapthisteq(I8, 'ClipLimit', CL);

V8 = fibermetric(J8, 'ObjectPolarity', 'bright');
V8_1 = V8>=0.2;


close all
%% Show images
subplot(3,1,1)
imshow(image,[min(I) max(I)])
title('Original 16-bit image')
subplot(3,1,2)
imshow(J,[min(min(J)) max(max(J))])
title('CLAHE Results for 16-bit image')
subplot(3,1,3)
imshow(V1)
title('Vesselness for 16-bit CLAHE Result')


figure,
subplot(3,1,1)
imshow(I8,[min(img8) max(img8)])
title('Original 8-bit image')
subplot(3,1,2)
imshow(J8,[min(min(J8)) max(max(J8))])
title('CLAHE Results for 8-bit image')
subplot(3,1,3)
imshow(V8)
title('Vesselness for 8-bit CLAHE Result')
