%% Locate images
close all; clear all; clc;
location = uigetdir;
newLocation = strcat(location, '\Demons Registered\');
mkdir(newLocation)

    %% Register Images
    fixed = imread([location '\0840.tif']); %reference image
%     fixedGPU = gpuArray(fixed); %Create gpuArray 
%     fixed1 = getMatrixOutliers(fixed);
%     fixed2 = fixed(find(fixed>0));
%     high = max(fixed2);
%     low = min(fixed2);
    
    cd(newLocation);
    
%     figure, subplot(4,4,1);
    for i = 120:120:1680
        imgpath = [location '\' sprintf('%04d.tif',i)];
        moving = imread(imgpath);
%         movingGPU = gpuArray(moving); %Create gpuArray


        %Perform Registration
%         [D, movingReg] = imregdemons(movingGPU, fixedGPU,200);
%         [D, movingReg] = imregdemons(movingGPU, fixedGPU);
        [~,movingReg] = imregdemons(moving,fixed,[500 400 200],'AccumulatedFieldSmoothing',3);
%         [~,movingReg] = imregdemons(movingGPU,fixedGPU,[500 400 200],'PyramidLevels',3);
        
        %Bring registered image back to CPU
% %         registeredImage = gather(movingReg);
        
        imwrite(movingReg,sprintf('%04d.tif',i))
        
        %Display registered images
%         j = (i/120)+1;
%         subplot(4,4,j);
%         imshowpair(fixed, registeredImage,'blend');
    end
    
    %% 

%     figure, subplot(4,4,1);
%     %Show all 15 newly registeted images on subplot 
%     for k = 0000:120:1680
%         path = [newLocation '\' sprintf('%04d.tif',k)];
%         image = imread(path);
%         imageGPU = gpuArray(image);
%         m = (k/120)+1;
%         subplot(4,4,m);
%         imshow(imageGPU, [low high]);
%     end