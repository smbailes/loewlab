path = uigetdir;
location = strcat(path, '\');


%strt = 1; %patient number         
    I_mat{i} = imread([location sprintf('%d.tif',strt)]);    % Read each image into I_mat
    newLocation = strcat(location, '\', 'Cropped Rect (individual)');
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
    imwrite(newCrop, sprintf('P%d_Cropped.tif',strt)); 
    
    strt=strt+1;            % Go to next image (for cropped), HAS TO BE CHANGED TO INCREMENT BY 1
close all

% newLocation = strcat(location, '\', 'Cropped Rect (individual)');
% mkdir(newLocation);
% 
% %% Apply Crop to All Registered Images based on User Drawn Input
% % cd ([homedir 'Registered/' ptID 'Registered/']);
% k = 1;
% for n = 0:120:1680
%     image = I_mat{k};
%     fprintf('%d\n',k);
%     I = getMatrixOutliers(image);
%     nonzero = I(find(I>0));
%     h = max(nonzero);
%     l = min(nonzero);
% 
%     figure
%     set(gcf,'units','inches', 'Position',[4 2 10 8])
%     imshow(image,[l h]);
% 
%     rect = imrect();
%     xy = wait(rect);
%     
%     newCrop = imcrop(image, xy);
%     
%     cd(newLocation)    
%     imwrite(newCrop, sprintf('%04d.tif',n));  
%     k = k+1;
%     
% end
