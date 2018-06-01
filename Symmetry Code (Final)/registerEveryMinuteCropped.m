function registerEveryMinuteCropped(user, ptID) 

% Takes User Input for Patient Directory and applies a Rectangular crop to
% all images in the Patient Directory
%% Inputs
%Input dialog box 
if nargin == 0
    ptID = patientselect;% obtain first value in answer matrix
    user = userselect;          % Dialog Box for user selection
end 

a = strcmp(user,'Aidan');   % Compare user input to 'Aidan"
b = strcmp(user,'Lovelace');
if (a == 1 && b == 0)                 % If Aidan was selected
    location = (['/Users/AidanMurray/Data/GWBox/GRP_Loew-Doc/NadaKamona/Clinic Patients/Manual Crop/' ptID 'Cropped/Cropped Image']);
    newLocation = (['/Users/AidanMurray/Data/GWBox/GRP_Loew-Doc/NadaKamona/Clinic Patients/Registered/' ptID 'Registered/Cropped/']);
elseif (a == 0 && b == 0)
    location = (['\Users\shann\Box\GRP_Loew-Doc\NadaKamona\Clinic Patients\Manual Crop\' ptID 'Cropped\Cropped Image\']);
    newLocation = (['\Users\shann\Box\GRP_Loew-Doc\NadaKamona\Clinic Patients\Registered\' ptID 'Registered\Cropped\']);
elseif (a == 0 && b == 1)
    location = (['D:\BreastPatients\Manual Crop\' ptID 'Cropped\']);
    newLocation = (['D:\BreastPatients\Registered\' ptID 'Registered\']);
end     
      

mkdir(newLocation)
        
%% Image Crop: Each Image
cd(newLocation)

stay = imread([location '/0000.tif']);
imwrite(stay, '0000.tif');

for i = 0120:120:1680
    newFile = [location '/' sprintf('%04d.tif',i)];
    newImage = imread(newFile);
    figure
    set(gcf,'units','inches', 'Position',[4 2 10 8])
    imshow(newImage,[]);
%     imcontrast
    
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumStepLength = 0.05;
    optimizer.MaximumIterations = 100; %SETTING FOR NEW MATLAB

    registeredImage = imregister(newImage,stay,'similarity',optimizer,metric);
    imshow(registeredImage,[]);
    imwrite(registeredImage,sprintf('%04d.tif',i))
    close;
end