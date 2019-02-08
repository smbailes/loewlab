clear all;
location = uigetdir;
location = strcat(location, '\');
ptID = 'IRSTXXX';

t1 = imread([location sprintf('%04d.tif',0)]);
t15 = imread([location sprintf('%04d.tif',1680)]);

I = getMatrixOutliers(t1);
I1 = I(find(I>0));
h = max(I1);
l = min(I1);

%% Get stats on the images
range1 = max(max(t1)) - min(min(t1));
range15 = max(max(t15)) - min(min(t15));
stdev1 = std2(t1);
stdev15 = std2(t15);
avg1 = mean2(t1);
avg15 = mean2(t15);

%% 
% eq = adapthisteq(t1, 'NumTiles',[15 15],'NBins',1000,'Range','original');
% le = min(eq(find(eq>0)));
% he = max(eq(find(eq>0)));
% figure, imshow(eq, [le he]);

%% Subtract images
diff = t1 - t15;
nonzero = diff(find(diff>0));
high = max(nonzero);
low = min(nonzero);
figure,
imshow(diff, [low high])
diffRange = high-low;
diffAvg = mean2(diff)
diffStd = std2(diff)

%% Fibermetric 
V = fibermetric(diff,'ObjectPolarity','dark','StructureSensitivity',50);
figure,
imshow(V)

%{
%Extracts indices of fibers
[r c] = size(V)
k = 1;
for i = 1:r
    for j = 1:c
        if(V(i, j)>0.25 && V(i,j)<0.85)
            ind{k,:} = [i j];
            k = k+1;
            t1(i,j) = 0;
        else
            V(i,j) = 0;
        end
    end
end
figure, imshow(V)
figure, imshow(t1, [l h])

%% Fibermetric on 1
V1 = fibermetric(t1, 10,'ObjectPolarity','bright','StructureSensitivity',12);
figure, 
imshow(V1)

[r c] = size(V1)
k = 1;
for i = 1:r
    for j = 1:c
        if(V1(i, j)>0.25 && V(i,j)<0.85)
            ind1{k,:} = [i j];
            k = k+1;
        end
    end
end


%% Fibermetric on 15
V15 = fibermetric(t15, 10,'ObjectPolarity','bright','StructureSensitivity',12);
figure, 
imshow(V15)

[r c] = size(V15)
k = 1;
for i = 1:r
    for j = 1:c
        if(V15(i, j)>0)
            ind15{k,:} = [i j];
            k = k+1;
        end
    end
end
%% Set results of fibermetric to zero

%% Try DBSCAN on new images? 
%% DBSCAN Parameters
minPts = 10; percent = 15; 
s = 1;  
epsilon = 6.5;

%% DBSCAN
%[ClustStruct, ClustData] = symmetry_cluster1(diff, epsilon, minPts, ptID, s, percent);

%}

