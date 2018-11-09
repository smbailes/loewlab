newLocation = uigetdir;
figure
%Find max and min for contrast 
path = [newLocation '/' sprintf('0000 - P11.tif')];
image = imread(path);
% I = getMatrixOutliers(image);
% I_nonzero = I(find(I>0));
% max = max(I_nonzero);
% min = min(I_nonzero);
    imshow(image, []);
