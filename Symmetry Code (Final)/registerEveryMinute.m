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

%% Register Images

stay = imread([location '\0840.tif']);
imwrite(stay,'0840.tif');
cd(newLocation);

for i = 0000:120:1680
    newFile = [location '\' sprintf('%04d.tif',i)];
    newImage = imread(newFile);
    nonzero = newImage(find(newImage>0));
%     figure
%     set(gcf,'units','inches', 'Position',[4 2 10 8])
%     imshow(newImage,[]);
%     imcontrast
    
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumStepLength = 0.05;
    optimizer.MaximumIterations = 100; %SETTING FOR NEW MATLAB

    registeredImage = imregister(newImage,stay,'affine',optimizer,metric);
%     figure
%     set(gcf,'units','inches', 'Position',[4 2 10 8])
%     imshow(registeredImage, []);
    imwrite(registeredImage,sprintf('%04d.tif',i))
    
end
figure, set(gcf,'units','inches', 'Position',[4 2 10 8]), subplot(4,4,1); 
for j = 1:15
    for k = 0000:120:1680
        path = [newLocation '\' sprintf('%04d.tif',k)];
        image = imread(path);
    end
   subplot(4,4,j);
   imshow(image, [min(nonzero) max(nonzero)]);
end 
