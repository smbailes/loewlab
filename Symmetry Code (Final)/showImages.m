newLocation = uigetdir;
figure
%Find max and min for contrast 
path = [newLocation '\' sprintf('%04d.tif',0)];
image = imread(path);
% I = getMatrixOutliers(image);
% I_nonzero = I(find(I>0));
% max = max(I_nonzero);
% min = min(I_nonzero);

for k = 0000:120:1680
    path = [newLocation '\' sprintf('%04d.tif',k)];
    image = imread(path);
    i = (k/120)+1;
    subplot(4,4,i);
    imshow(image, []);
end