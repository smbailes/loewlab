%% Loading Fixed Image
ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
I = imread([dir '/' ptID]);

fixed = I;

JID = input('Enter image name for the moving image to open: ','s');
JID = strcat(JID,'.tif');
dir = uigetdir;
newImage = imread([dir '/' JID]);

%% Resize Images

fprintf('Pick right nipple. \n');
figure, imshow(newImage, []), title('Right Nipple Visible')
[X1,Y1] = ginput(1);

fprintf('Pick right nipple. \n');
figure, imshow(I, []), title('Right Nipple Reference')
[X2,Y2] = ginput(1);

val = mean(mean(newImage));

% Cropping

hshift = round(X2 - X1);
if hshift>0
    newImage = padarray(newImage, [0 abs(hshift)], val, 'pre');
    newImage = imcrop(newImage, [0 0 640 512]);
else
    newImage = padarray(newImage, [0 abs(hshift)], val, 'post');
    [r c] = size(newImage);
    newImage = imcrop(newImage, [(c-640) (512-r) (c) (r)]);
end

vshift = round(Y2 - Y1);
if vshift>0
    newImage = padarray(newImage, abs(vshift), val, 'pre');
    newImage = imcrop(newImage, [0 0 640 512]);
else
    newImage = padarray(newImage, abs(vshift), val, 'post');
    [r c] = size(newImage);
    newImage = imcrop(newImage, [(640-c) (r-511) c r]); 
end

figure(1), imshow(newImage)
figure(2), imshow(I, [])


%% Affine Registration

[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.01;
optimizer.MaximumIterations = 100; %SETTING FOR NEW MATLAB
%optimizer.Epsilon = 1.5e-8; %epsilon = max diff between 2 imgs
optimizer.GrowthFactor = 1.05;

registeredImage = imregister(newImage,fixed,'affine',optimizer,metric);
imwrite(registeredImage,sprintf('VisibleAff.tif',i))

fprintf('Finished Affine Registration\n');

%% Demons Registration
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

%% Show Affine Registration

affID = input('Enter image name for the affine registered image to open: ','s');
affID = strcat(affID,'.tif');
dir = uigetdir;
affI = imread([dir '/' affID]);

affI = imshow(affI, []);
%go on to manual segmentation