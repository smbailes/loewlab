function registerEveryMinute(user, ptID) 

% Takes User Input for Patient Directory and applies a Rectangular crop to
% all images in the Patient Directory
%% Inputs
%Input dialog box 
%     [location, ptID] = pathfinder;
%     ptID = patientselect;
    location = uigetdir;
    newLocation = strcat(location, '\Registered\');
    mkdir(newLocation)

%% Register Images and Show Alignment 
ref = imread([location '\0840.tif']); %reference image
ref1 = getMatrixOutliers(ref);
ref_nonzero = ref1(find(ref1>0));
high = max(ref_nonzero);
low = min(ref_nonzero);

cd(newLocation);
figure, subplot(4,4,1);
for i = 0000:120:1680
    newFile = [location '\' sprintf('%04d.tif',i)];
    I = imread(newFile);
    I_outliers = getMatrixOutliers(I);
%     figure
%     set(gcf,'units','inches', 'Position',[4 2 10 8])
%     imshow(newImage,[]);
%     imcontrast
%     figure

    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumStepLength = 0.05;
    optimizer.MaximumIterations = 150; %SETTING FOR NEW MATLAB

    registeredImage = imregister(I_outliers,ref,'affine',optimizer,metric);
    imwrite(registeredImage,sprintf('%04d.tif',i))
    j = (i/120)+1;
    subplot(4,4,j);
    imshowpair(ref, registeredImage, 'Scaling', 'joint');
end
  
figure, subplot(4,4,1);
%Show all 15 newly registeted images on subplot 
for k = 0000:120:1680
    path = [newLocation '\' sprintf('%04d.tif',k)];
    image = imread(path);
    m = (k/120)+1;
    subplot(4,4,m);
    imshow(image, [low high]);
end
end