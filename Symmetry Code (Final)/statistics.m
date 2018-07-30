clear all, clc
%% Patient Selection
    [location, ptID] = pathfinder; 

%% Image Input
    
    % Read 14 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end

%% Get statistics 

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

