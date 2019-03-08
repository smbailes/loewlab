ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
I = imread([dir '\' ptID]);

fixed = I_mat{8};
cd(newLocation);

%Affines Registration
JID = input('Enter image name for the moving image to open: ','s');
JID = strcat(JID,'.tif');
dir = uigetdir;
newImage = imread([dir '\' JID]);
    
[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumStepLength = 0.04;
optimizer.MaximumIterations = 100; %SETTING FOR NEW MATLAB

registeredImage = imregister(newImage,fixed,'affine',optimizer,metric);
imwrite(registeredImage,sprintf('VisibleAff.tif',i))

fprintf('Finished Affine Registration\n');

%Demons Registration
fixedGPU = gpuArray(fixed); %Create gpuArray
movID = input('Enter image name for the moving image to open: ','s');
movID = strcat(movID,'.tif');
dir = uigetdir;
moving = imread([dir '\' movID]);
%         movingGPU = gpuArray(moving); %Create gpuArray
[~,movingReg] = imregdemons(moving,fixed,[500 400 200],'AccumulatedFieldSmoothing',3);
%         [~,movingReg] = imregdemons(movingGPU,fixedGPU,[500 400 200],'AccumulatedFieldSmoothing',3);
        
        %Bring registered image back to CPU
%         registeredImage = gather(movingReg);
        
imwrite(movingReg,sprintf('VisibleDem.tif',i));
fprintf('Finished Demons Registration\n');

pause 
close all;
