clear all, close all
%% Patient Selection
    [location, ptID] = pathfinder; 

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end
    
%% Find statistics 

[r c] = size(I_mat{1}); 

for i = 1:15
    I1 = I_mat{i}(find(I_mat{i}>0));
    highcol = max(I1);
    high(i) = max(highcol); 
    lowcol = min(I1);
    low(i) = min(lowcol);
    range(i) = high(i) - low(i);
    average(i) = mean2(I1);
    stdev(i) = std2(I1);
end 

BreastInfo(1:1) = struct('Volunteer',ptID,'Average',average,'StandardDeviation',stdev,'Range',range,'HighValue',high,'LowValue',low)

%% Find stats for each side
image = I_mat{8};
I = image(find(image>0));
h = max(I);
l = min(I);
figure, imshow(I_mat{8},[l h]);
[x,y] = ginput(1);




