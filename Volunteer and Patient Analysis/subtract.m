clear all;
location = uigetdir;
location = strcat(location, '\');
ptID = 'IRSTXXX';

t1 = imread([location sprintf('%04d.tif',0)]);
t15 = imread([location sprintf('%04d.tif',1680)]);

%% Get stats on the images
range1 = max(max(t1)) - min(min(t1));
range15 = max(max(t15)) - min(min(t15));
stdev1 = std2(t1);
stdev15 = std2(t15);
avg1 = mean2(t1);
avg15 = mean2(t15);

%% 
eq = adapthisteq(t1, 'NumTiles',[15 15],'NBins',1000,'Range','original');
le = min(eq(find(eq>0)));
he = max(eq(find(eq>0)));
figure, imshow(eq, [le he]);

%% Subtract images
diff = t1 - t15;
nonzero = diff(find(diff>0));
high = max(nonzero);
low = min(nonzero);
figure,
imshow(diff, [low high])
diffRange = high-low;
diffAvg = mean2(diff);
%% Fibermetric 
V = fibermetric(diff,10,'ObjectPolarity','bright','StructureSensitivity',12);
figure,
imshow(V)

%% Try DBSCAN on new images? 
%% DBSCAN Parameters
minPts = 10; percent = 15; 
s = 1;  
epsilon = 6.5;

%% DBSCAN
%[ClustStruct, ClustData] = symmetry_cluster1(diff, epsilon, minPts, ptID, s, percent);




