fixed = I_mat{8};
cd(newLocation);

if strcmp(in, 'a')%Affine registration
    for i = 0:120:1680
    newFile = [location '/' sprintf('%04d.tif',i)];
    newImage = imread(newFile);
    
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumStepLength = 0.04;
    optimizer.MaximumIterations = 100; %SETTING FOR NEW MATLAB

    registeredImage = imregister(newImage,fixed,'affine',optimizer,metric);
    imwrite(registeredImage,sprintf('%04d.tif',i))
    end
    fprintf('Finished Affine Registration\n');
elseif strcmp(in, 'd') %demons registration
    
    fixedGPU = gpuArray(fixed); %Create gpuArray
    for i = 0:120:1680
        imgpath = [location '\' sprintf('%04d.tif',i)];
        moving = imread(imgpath);
%         movingGPU = gpuArray(moving); %Create gpuArray
        [~,movingReg] = imregdemons(moving,fixed,[500 400 200],'AccumulatedFieldSmoothing',3);
%         [~,movingReg] = imregdemons(movingGPU,fixedGPU,[500 400 200],'AccumulatedFieldSmoothing',3);
        
        %Bring registered image back to CPU
%         registeredImage = gather(movingReg);
        
        imwrite(movingReg,sprintf('%04d.tif',i))
    end 
    fprintf('Finished Demons Registration\n'); 
end

%% Show newly registered images
n = 0;
for j = 1:15
    I_reg{j} = imread([newLocation sprintf('%04d.tif',n)]);
    n = n+120;
end

for p = 1:15
    figure
    imshow(I_mat{p}, [l h]);
end
pause 
close all;

%% Crop Images
RegisteredCropScript;
pause, close all