%% Setup: Find directory, open image
close all; clear all; clc;
location = uigetdir;
location = strcat(location, '\');

%% Read in Images
strt = 0;
for i=1:15          
    I_mat{i} = imread([location sprintf('%04d.tif',strt)]);    % Read each image into I_mat
    strt=strt+120;            % Go to next image (for cropped), HAS TO BE CHANGED TO INCREMENT BY 1
end
newLocation = strcat(location, '\Registered\');
mkdir(newLocation)

%% Check if small/large breasts
image = I_mat{8};
I = getMatrixOutliers(image);
nonzero = I(find(I>0));
h = max(nonzero);
l = min(nonzero);

figure
set(gcf,'units','inches', 'Position',[4 2 10 8])
imshow(image,[l h]);
pause(2);
close;

in = input('Is the breast small or large? Enter s/l: ','s');

%% Register Images
fixed = I_mat{8};
cd(newLocation);

if strcmp(in, 's')%Affine registration
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
elseif strcmp(in, 'l') %demons registration
    fixedGPU = gpuArray(fixed); %Create gpuArray
    cd(newLocation);

    for i = 0:120:1680
        imgpath = [location '\' sprintf('%04d.tif',i)];
        moving = imread(imgpath);
        movingGPU = gpuArray(moving); %Create gpuArray

        [~,movingReg] = imregdemons(movingGPU,fixedGPU,[500 400 200],'AccumulatedFieldSmoothing',3);
        
        %Bring registered image back to CPU
        registeredImage = gather(movingReg);
        
        imwrite(registeredImage,sprintf('%04d.tif',i))
    end 
    fprintf('Finished Demons Registration\n'); 
end
