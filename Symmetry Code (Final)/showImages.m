newLocation = uigetdir;
%Find max and min for contrast 
path = [newLocation '\' sprintf('%04d.tif',840)];
image = imread(path);
I = getMatrixOutliers(image);
I_nonzero = I(find(I>0));
h = max(I_nonzero);
l = min(I_nonzero);

figure
for k = 120:120:1680
    path = [newLocation '\' sprintf('%04d.tif',k)];
    image = imread(path);
    i = (k/120)+1;
    subplot(4,4,i);
    imshow(image, [l h]);
end