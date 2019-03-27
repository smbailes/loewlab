close all
ptID = 'IRST019';
i = 0;
newLocation = uigetdir;
%Find max and min for contrast 
path = [newLocation '\' sprintf('%s-%04d.tif',ptID,i)];
image = imread(path);
% I = getMatrixOutliers(image);
I_nonzero = image(find(image>0));
h = max(I_nonzero);
l = min(I_nonzero);

% figure
%Show in subplots
% for k = 120:120:1680
%     path = [newLocation '\' sprintf('%04d.tif',k)];
%     image = imread(path);
%     i = (k/120)+1;
%     subplot(4,4,i);
%     imshow(image, [l h]);
% end

for k = 0:120:1680
    path = [newLocation '\' sprintf('%04d.tif',k)];
    image = imread(path);
    figure
    imshow(image,[l h]);
end
