%% Setup: Find directory, open image
% Dependencies
% - getMatrixOutliers 
% - RegisteredCropScript 
% - showImages 

close all; clear all; clc;
location = uigetdir;
location = strcat(location, '\');

%% Read in Images
strt = 0;
for i=1:15          
    I_mat{i} = imread([location sprintf('%04d.tif',strt)]);    % Read each image into I_mat
    strt=strt+120;            % Go to next image (for cropped), HAS TO BE CHANGED TO INCREMENT BY 1
end
newLocation = strcat(location, 'Registered 2\');
mkdir(newLocation)

%% Show images to determine whether to use affine or demons registration
image = I_mat{8};
I = getMatrixOutliers(image);
nonzero = I(find(I>0));
h = max(nonzero);
l = min(nonzero);

for i = 1:15
    figure
    imshow(I_mat{i}, [l h]);
end
pause
close all
%% Register Images
% Use demons for patients that move more or with larger breasts

in = input('Affine or demons? Enter a/d: ','s');
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
    % Check to see if computer has GPU. If it does, use GPU. 
    g = input('Does this computer have a GPU? (Y/N):', 's'); 
    if strcmp(g, 'Y')
        fixedGPU = gpuArray(fixed); %Create gpuArray
        for i = 0:120:1680
            imgpath = [location '\' sprintf('%04d.tif',i)];
            moving = imread(imgpath);
            movingGPU = gpuArray(moving); %Create gpuArray

            [~,movingReg] = imregdemons(movingGPU,fixedGPU,[500 400 200],'AccumulatedFieldSmoothing',3);

            %Bring registered image back to CPU
            registeredImage = gather(movingReg);

            imwrite(registeredImage,sprintf('%04d.tif',i))
        end 
    else
         for i = 0:120:1680
            imgpath = [location '\' sprintf('%04d.tif',i)];
            moving = imread(imgpath);
            [~,registeredImage] = imregdemons(moving,fixed,[500 400 200],'AccumulatedFieldSmoothing',3);
            imwrite(registeredImage,sprintf('%04d.tif',i))
         end 
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
    imshow(I_reg{p}, [l h]);
end
pause 
close all;

%% Crop Images
RegisteredCropScript;
pause, close all