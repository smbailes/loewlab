path = uigetdir;
location = strcat(path, '\');


strt = 0;
for i=1:15          
    I_mat{i} = imread([location sprintf('%04d.tif',strt)]);    % Read each image into I_mat
    strt=strt+120;            % Go to next image (for cropped), HAS TO BE CHANGED TO INCREMENT BY 1
end

newLocation = strcat(location, '\', 'Cropped Rect');
mkdir(newLocation);

%% Apply Crop to All Registered Images based on User Drawn Input
% cd ([homedir 'Registered/' ptID 'Registered/']);
image = I_mat{8};
I = getMatrixOutliers(image);
nonzero = I(find(I>0));
h = max(nonzero);
l = min(nonzero);

figure
set(gcf,'units','inches', 'Position',[4 2 10 8])
imshow(image,[l h]);

rect = imrect();
xy = wait(rect);


% close;
n = 1;
for k = 0:120:1680    
    newCrop = imcrop(I_mat{n}, xy);
    
    cd(newLocation)    
    imwrite(newCrop, sprintf('%04d.tif',k));  
    n = n+1;
end
close;
