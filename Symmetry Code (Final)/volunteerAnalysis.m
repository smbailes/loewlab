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

%% Split into right and left breast 
image = I_mat{8};
I = image(find(image>0));
h = max(I);
l = min(I);
figure, imshow(I_mat{8},[l h]);
title('Select Patients Right Breast')

rect = imrect();
binaryImage = rect.createMask();
rightBW = uint16(binaryImage);
leftBW = 1-rightBW;

for i = 1:15
    rightBreast{i} = I_mat{i}.*rightBW;
    leftBreast{i} = I_mat{i}.*leftBW;
end 

%% Find statistics for right and left

for i = 1:15
    I1 = rightBreast{i}(find(rightBreast{i}>0));
    highcolRight = max(I1);
    highRight(i) = max(highcolRight); 
    lowcolRight = min(I1);
    lowRight(i) = min(lowcolRight);
    rangeRight(i) = highRight(i) - lowRight(i);
    averageRight(i) = mean2(I1);
    stdevRight(i) = std2(I1);
end 


for i = 1:15
    I1 = leftBreast{i}(find(leftBreast{i}>0));
    highcolLeft = max(I1);
    highLeft(i) = max(highcolLeft); 
    lowcolLeft = min(I1);
    lowLeft(i) = min(lowcolLeft);
    rangeLeft(i) = highLeft(i) - lowLeft(i);
    averageLeft(i) = mean2(I1);
    stdevLeft(i) = std2(I1);
end 

BreastRightInfo(1:1) = struct('Side','Right','Average',averageRight,'StandardDeviation',stdevRight,'Range',rangeRight,'HighValue',highRight,'LowValue',lowRight)
BreastLeftInfo(1:1) = struct('Side','Left','Average',averageLeft,'StandardDeviation',stdevLeft,'Range',rangeLeft,'HighValue',highLeft,'LowValue',lowLeft)


%% Look at histogram across vessel
figure, subplot(2,1,1), imshow(I_mat{8},[l h]);
title('Select rectangle going across a vessel')
rect = imrect();
binaryImage = rect.createMask();
acrossVessel = I_mat{8}.*(uint16(binaryImage));
nonzero = acrossVessel(find(acrossVessel>0));
subplot(2,1,2), plot(nonzero,'-o')

again = questdlg('Try again?', 'Yes', 'No');

while strcmp(again,'Yes') == 1
    
    figure, subplot(2,1,1), imshow(I_mat{8},[l h]);
    title('Select rectangle going across a vessel')
    rect = imrect();
    binaryImage = rect.createMask();
    acrossVessel = I_mat{8}.*(uint16(binaryImage));
    nonzero = acrossVessel(find(acrossVessel>0));
    subplot(2,1,2), plot(nonzero,'-o')

    again = questdlg('Try again?','Yes','No');
end 
