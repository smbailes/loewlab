%function registerEveryMinute(user, ptID) 

    % Takes User Input for Patient Directory and applies a Rectangular crop to
    % all images in the Patient Directory
    %% Inputs
    %Input dialog box 
    %     [location, ptID] = pathfinder;
    %     ptID = patientselect;
         location = uigetdir;
         ptID = input('Enter patient number: ','s'); %Request patient image name
    %     ptID = strcat(ptID,'.tif'); 
        newLocation = strcat(location, '\Registered\');    
        location = strcat(location,'\', ptID); 
        
         mkdir(newLocation)
         

    %% Register Images and Show Alignment 
    %path1 = uigetfile('.tif');
    %ref = imread(path1);

    ref = imread([location '-0015.tif']); %reference image
    ref1 = getMatrixOutliers(ref);
    ref_nonzero = ref1(find(ref1>0));
    high = max(ref_nonzero);
    low = min(ref_nonzero);

    cd(newLocation);
    figure, subplot(4,4,1);
    for i = 0001:1:0030
        newFile = [location '-' sprintf('%04d.tif',i)];
        I = imread(newFile);
%         I_outliers = getMatrixOutliers(I);
%         figure
%         set(gcf,'units','inches', 'Position',[4 2 10 8])
%         imshow(newImage,[]);
%         imcontrast
%         figure

        [optimizer, metric] = imregconfig('monomodal');
        optimizer.MaximumStepLength = 0.05;
        optimizer.MaximumIterations = 150; %SETTING FOR NEW MATLAB

        registeredImage = imregister(I,ref,'affine',optimizer,metric);
        imwrite(registeredImage,sprintf('%04d.tif',i))
        j = (i/120)+1;
        subplot(4,4,j);
        imshowpair(ref, registeredImage, 'Scaling', 'joint');
    end
    %% 

    figure, subplot(4,4,1);
    %Show all 15 newly registeted images on subplot 
    for k = 0001:1:0030
        path = [newLocation, sprintf('%04d.tif',k)];
        image = imread(path);
        m = (k/120)+1;
        subplot(4,4,m);
        imshow(image, [low high]);
    end

% %% Registration of JPEG
% 
%     ref1 = getMatrixOutliers(ref);
%     ref_nonzero = ref1(find(ref1>0));
%     high = max(ref_nonzero);
%     low = min(ref_nonzero);
% 
%     cd(newLocation);
%     figure, subplot(4,4,1);
%     for i = 1000:1000:4000
%         newFile = [location '-' sprintf('%04d.tif',i)];
%         I = imread(newFile);
% %         I_outliers = getMatrixOutliers(I);
% %         figure
% %         set(gcf,'units','inches', 'Position',[4 2 10 8])
% %         imshow(newImage,[]);
% %         imcontrast
% %         figure
% 
%         [optimizer, metric] = imregconfig('monomodal');
%         optimizer.MaximumStepLength = 0.05;
%         optimizer.MaximumIterations = 150; %SETTING FOR NEW MATLAB
% 
%         registeredImage = imregister(I,ref,'affine',optimizer,metric);
%         imwrite(registeredImage,sprintf('%04d.tif',i))
%         j = (i/120)+1;
%         subplot(4,4,j);
%         imshowpair(ref, registeredImage, 'Scaling', 'joint');
%     end
    
%% Registration of JPEG 2

vis = input('Enter image name you want to open: ','s'); 
vis = strcat(vis,'.tif'); 

dir1 = uigetdir; 
Ivis = imread([dir1 '\' vis]); 

ref = input('Enter image name you want to open: ','s'); 
ref = strcat(ref,'.tif'); 

dir2 = uigetdir; 
Iref = imread([dir2 '\' ref]); 

fprintf('Pick right nipple. \n');
figure, imshow(Ivis, []), title('Right Nipple Visible')
[X1,Y1] = ginput(1);

fprintf('Pick right nipple. \n');
figure, imshow(Iref, []), title('Right Nipple Reference')
[X2,Y2] = ginput(1);

val = mean(mean(Ivis));

% Cropping

hshift = round(X2 - X1);
if hshift>0
    Ivis = padarray(Ivis, [0 hshift], val, 'pre');
    Ivis = imcrop(Ivis, [0 0 640 512]);
else
    Ivis = padarray(Ivis, [0 abs(hshift)], val, 'post');
    [r c] = size(Ivis);
    Ivis = imcrop(Ivis, [r c 512 640]);
end

%% 

vshift = round(Y2 - Y1);
if vshift>0
    Ivis = padarray(Ivis, vshift, val, 'pre');
    Ivis = imcrop(Ivis, [0 0 512 640]);
else
    Ivis = padarray(Ivis, abs(vshift), val, 'post');
    [r c] = size(Ivis);
    Ivis = imcrop(Ivis, [0 c 512 640]);
end

figure(1), imshow(Ivis, [])
figure(2), imshow(Iref, [])
%% Scaling Down
% 
% fprintf('Pick right nipple. \n');
% figure, imshow(Ivis, []), title('Right Nipple Visible')
% [Xr1,Yr1] = ginput(1);
% 
% fprintf('Pick right nipple. \n');
% figure, imshow(Iref, []), title('Right Nipple Reference')
% [Xr2,Yr2] = ginput(1);
% 
% fprintf('Pick left nipple. \n');
% figure, imshow(Ivis, []), title('Left Nipple Visible')
% [Xl1,Yl1] = ginput(1);
% 
% fprintf('Pick left nipple. \n');
% figure, imshow(Iref, []), title('Left Nipple Reference')
% [Xl2,Yl2] = ginput(1);
%  

% visdist = sqrt((Xr1^2)-(Xl1^2));
% refdist = sqrt((Xr2^2)-(Xl2^2));
% 
% if visdist > refdist
%     ratio = refdist/visdist;
%     Ivis = imresize(Ivis, ratio);
% else 
%     ratio = refdist/visdist;
%     Ivis = imresize(Ivis, ratio);
%     Ivis = imcrop(Ivis, [0 0 640 512]);
% end 


%% Compare

II = imshowpair(Ivis,Iref)
%% 
% 
% hold on
% % Save the handle for later use 
% h = imshow(Ivis); 
% hold off
% 
% [M,N] = size(Ivis); 
% block_size = 100; 
% P = ceil(M / block_size); 
% Q = ceil(N / block_size); 
% alpha = checkerboard(block_size, P, Q) > 0; 
% alpha = alpha(1:M, 1:N); 
% set(h, 'AlphaData', alpha);


