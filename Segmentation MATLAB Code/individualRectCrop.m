path = uigetdir;
location = strcat(path, '/');


        
strt=12
i=strt
    I_mat{i} = imread([location sprintf('0000 - P11',strt)]);    % Read each image into I_mat %,strt
    
    newLocation = strcat(location, '/', '0000 - P11');
    mkdir(newLocation);
    
    image = I_mat{i};
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
    imwrite(newCrop, sprintf('Pat.tif',strt));  
    strt=strt+1;            % Go to next image (for cropped), HAS TO BE CHANGED TO INCREMENT BY 1
close all
