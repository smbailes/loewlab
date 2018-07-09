path = uigetdir;
location = strcat(path, '\');


strt = 0;
for i=1:15          
    I_mat{i} = imread([location sprintf('%04d.tif',strt)]);    % Read each image into I_mat
    strt=strt+120;            % Go to next image (for cropped), HAS TO BE CHANGED TO INCREMENT BY 1
end

newLocation = strcat(location, '\', 'Cropped Rect (individual)');
mkdir(newLocation);

%% Apply Crop to All Registered Images based on User Drawn Input
% cd ([homedir 'Registered/' ptID 'Registered/']);
k = 1;
for n = 0:120:1680
    image = I_mat{k};
    fprintf('%d\n',k);
    I = getMatrixOutliers(image);
    nonzero = I(find(I>0));
    h = max(nonzero);
    l = min(nonzero);

    figure
    set(gcf,'units','inches', 'Position',[4 2 10 8])
    imshow(image,[l h]);

    rect = imrect();
    xy = wait(rect);
    
    newCrop = imcrop(image, xy);
    
    cd(newLocation)    
    imwrite(newCrop, sprintf('%04d.tif',n));  
    k = k+1;
    
end
