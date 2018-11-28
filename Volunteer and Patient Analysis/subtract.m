clear all;
location = uigetdir;
location = strcat(location, '\');

t1 = imread([location sprintf('%04d.tif',0)]);
t15 = imread([location sprintf('%04d.tif',1680)]);
%% Subtract images
diff = t1 - t15;
nonzero = diff(find(diff>0));
high = max(nonzero);
low = min(nonzero);
figure,
imshow(diff, [low high])
high-low 
mean2(diff)
%% Fibermetric 
V = fibermetric(diff,10,'ObjectPolarity','bright','StructureSensitivity',12);
figure,
imshow(V)



