% Part 1 of the Code for Breast Segmentation
% Last Modified: Kate Fergusson
% 06/21/2017 4:32 PM

clear 
close all
clc

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
I = imread(ptID); %open the image, keeping it in 16-bits

figure, imshow(I,[]) %to help decide if it should be cropped or not
title('Image 1: Original Image')

cropopt=input('Does a crop need to occur?? y/n: ','s');
if cropopt == 'y'
    I=imcrop(I,[]); %cropping, if necessary
end

close Figure 1

tic
ed1 = edge(I,'canny');%finds canny edge
CCANNY = bwconncomp(ed1); %CCANNY is used later to find the boundaries of the canny edges

pxlength = 10; % lines to be removed if under 10 pixels long
boundededges = bwareaopen(ed1, pxlength); % act of removing said lines
smaller = bwmorph(boundededges,'skel',Inf); %move the canny edges down to skeleton
figure, imshow(smaller)
title('Image 2: Skeleton Canny Image, to find ellipses from')


figure;
imshow(I,[]); %produce an image to overlay the ellipses onto
title('Image with Ellipses')

%RIGHT SIDE
% override some default parameters
paramsr.minMajorAxis = 350;
paramsr.maxMajorAxis = 700;
paramsr.numBest = 12; %draws 12 ellipses
paramsr.rotation = 45; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
paramsr.rotationSpan = 35;
%paramsr.randomize = 0; %randomization component that may reduce changing of
%ellipses

% note that the edge (or gradient) image is used
bestFitsr = ellipseDetection(smaller, paramsr);
fprintf('Output %d best fits.\n', size(bestFitsr,1));


%ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 

%takes the information that was found of the ellipses and draws them;also
%keeping the information for each ellipse in a cell in qr(and later ql for those):
qr{1,length(bestFitsr)}=0;
for n=1:length(bestFitsr)
    qr{n} = ellipse(bestFitsr(n,3),bestFitsr(n,4),bestFitsr(n,5)*pi/180,bestFitsr(n,1),bestFitsr(n,2),'k');
end
%overriding parameters:
paramsl.minMajorAxis = 350;
paramsl.maxMajorAxis = 700;
paramsl.numBest = 12; %draws 12 ellipses
paramsl.rotation = 135; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
paramsl.rotationSpan = 35;
%paramsl.randomize = 0; %randomization component that may reduce changing of
%ellipses


%LEFT SIDE
bestFitsl = ellipseDetection(smaller, paramsl);
fprintf('Output %d best fits.\n', size(bestFitsl,1));

%ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
ql{1,length(bestFitsl)}=0;
for n=1:length(bestFitsl)
    ql{n} = ellipse(bestFitsl(n,3),bestFitsl(n,4),bestFitsl(n,5)*pi/180,bestFitsl(n,1),bestFitsl(n,2),'k');
end

%clear bestFitsl bestFitsr


%ellipse(semimajor axis of ra,semimajor axis of radius rb, semimajor axis
%of ang, center point x0, center pointy0, 'color')
%
% **The length of ra, rb, and ang should be the same.**

%the following is for finding the Hough circles, and then overlaying them
rmin = 75;
rmax = 125;
[centers, radii] = imfindcircles(smaller,[rmin rmax],'Sensitivity',0.97); 
%A higher 'Sensitivity' value sets the detection threshold lower and leads 
%to detecting more circles.
circles = viscircles(centers,radii,'LineWidth',1); 

if length(centers)>2   % keeps it to only 2 circles. (hopefully the correct ones)
    centers(3:end,:)=[];
end



% Part 2: find the boundaries of hot regions in newI.
% Taken from Nada's code

I_sobel = edge(I,'sobel');
se = strel('disk',2); %Create a Morphological structuring element, you change the shape used and diameter
I_sobel= imclose(I_sobel,se); % Morphologically close image using the se element

CC = bwconncomp(I_sobel); %Find connected components

%Find the length of each conncomp in CC and store them in an array
for n = 1:CC.NumObjects
    ConnCompLengths(n)= length(CC.PixelIdxList{n});
end
clear n;




[max1, indx1] = max(ConnCompLengths); %find longest components and its index
ConnCompLengths(indx1) = 0; % delete that max, to find the second maximum
[max2, indx2] = max(ConnCompLengths); %find second largest components

[r1,tmp] = ind2sub(size(I),CC.PixelIdxList{indx2}); %find the row/col index
r2 = max(r1); % top row for the first component
r1 = min(r1); % find the minimum row in one component of sobel

tmp2=max(tmp);
tmp1=min(tmp);


newI((1:r1),:) = 0; %Top rows will be deleted, to limit the area of interest between r1 and r2


[r,c] = size(I);

boundaries_hotpixel = zeros(size(I)); %empty image

%limit search area between r1 and r2, and find the boundaries in newI using
%the two loops below. 
for n = r1%:r2
    if (sum(newI(n,:))~=0) %the sum per row
        loc = find(newI(n,:)>0); %find non-zero values in that row
        boundaries_hotpixel(n,loc(1)) = 2^16; %take first non-zero pixel in each row, and set it to 1
        boundaries_hotpixel(n,loc(end)) = 2^16; %take the last non-zero pixel in each row
    end   
end

for n = 1%:c
    if (sum(newI((r1),n))~=0) %the sum per col :r2
        loc = find(newI((r1),n)>0); %find non-zero values in that col:r2
        boundaries_hotpixel((loc(1)+r1-1),n) = 2^16; %take first non-zero pixel in each col, and set it to 1
    end   
end


[img_y, img_x] = size(I);
xe1 = zeros(length(ql),[]);
ye1 = zeros(length(ql),[]);
ellipse_matrix = zeros(img_y,img_x);        %create matrix of zeros as big as the image



%draw in left breast ellipses pixel by pixel
for a = 1:length(ql)                %a,b,n,m are just used as counters in the for loops - delete at end of section
    e1 = ql{a};
    for b = 1:length(e1.XData)              %get x and y data from cell array of ellipses and round so we can use them as indices
        xe1(a,b) = e1.XData(b);
        xe1(a,b) = round(xe1(a,b));
        ye1(a,b) = e1.YData(b);
        ye1(a,b) = round(ye1(a,b));
        if xe1(a,b)<0.5
            xe1(a,b)=1;
        end
        if ye1(a,b)<0.5
            ye1(a,b)=1;
        end
    end
    for d = 1:length(xe1)
        ellipse_matrix(ye1(a,d),xe1(a,d)) = 1;      %fill in 1's wherever there is a point in the ellipse
    end
    for n = 1:img_y
        for m = 1:img_x
            if boundaries_hotpixel(n,m) == 0
                ellipse_matrix(n,m)=0;              %THIS IS NOT WORKING COMPLETELY YET!! Some pixels are still there
                %break;                              %if hotpixel = 0, then we want to remove that pixel of ellipse if its there
            end
        end
    end
end

%draw in right breast pixel by pixel
for a = 1:length(qr)                %a,b,n,m are just used as counters in the for loops - delete at end of section
    e2 = qr{a};
    for b = 1:length(e2.XData)
        xe2(a,b) = e2.XData(b);
        xe2(a,b) = round(xe2(a,b));
        ye2(a,b) = e2.YData(b);
        ye2(a,b) = round(ye2(a,b));   
        if xe2(a,b)<0.5
            xe2(a,b)=1;
        end
        if ye2(a,b)<0.5
            ye2(a,b)=1;
        end
    end
    for d = 1:length(xe2)
        ellipse_matrix(ye2(a,d),xe2(a,d)) = 1;      %fill in 1's wherever there is a point in the ellipse
    end
    for n = 1:img_y
        for m = 1:img_x
            if boundaries_hotpixel(n,m) == 0
                ellipse_matrix(n,m)=0;              %if hotpixel = 0, then we want to remove that pixel of ellipse if its there
            end
        end
    end
end




%_________________________________________________________________________%
 % Hot Pixel Part
bins = 2^16; %insert image bits here
[N,binlocation] = imhist(I,bins); %each count will has its own bin

perc = input('What is your desired percentage? '); 
lowerBound = .009*perc*numel(I); % numel(I): the number of pixels in the image
upperBound = .011*perc*numel(I);
global imTitle;
imTitle = sprintf('Hot Pixel, highest 5%%');
total = 0; 
for j = 1:numel(N)  
    total = total+ N(bins-j); %Tells us where the 1 percent of pixels start from but some of the bins still have zero in them
    if (total >= lowerBound)  
       break
    end 
end

k = 0;
hotPix(2) = 0;
for i = bins-j:bins %Keep all values that say 256 as 256. Matlab does not have a bin called 0 so there is an extra bin.
    if (N(i)~=0) %If the bin does not equal zero then move on to the next and add it to the hotPix list 
     k = k+1; 
     hotPix(k) = i;  %This contains the 1 percent of pixels that do not contain zero
     
    end 
end

newI = zeros(size(I)); %new image and make the entire image black
for m = 1:k
    o = find(I==hotPix(m)); %This looks for the hotPix in the image
    newI(o) = bins; %Make any value of the bightest 1 percent white
end
figure, imshow(newI)
title(imTitle)
clear o
clear m

boundaries_canny = zeros(size(I));
loc = find(newI>0); %match with the entire area in the p% from above part
for x = 1:length(loc) %go index by index
    for y = 1:CCANNY.NumObjects %go through each connected component in CC
        if find(CCANNY.PixelIdxList{y}==loc(x)) %check if location x is found in any of the connected components
            %curly brackets access the contents of the the cell
            boundaries_canny(CCANNY.PixelIdxList{y}) = 2^16; %if yes, then keep the ENTIRE component and store it in boundaries2
        end
    end
end
clear x
clear y


%Combine HotPixel and Canny
boundaries_hotpxcanny = boundaries_canny;

I_sobel = edge(I,'sobel');
se = strel('disk',2); %Create a Morphological structuring element, you change the shape used and diameter
I_sobel= imclose(I_sobel,se); % Morphologically close image using the se element

CC = bwconncomp(I_sobel); %Find connected components

%Find the length of each conncomp in CC and store them in an array
ConnCompLengths(2)=0;
for n = 1:CC.NumObjects
    ConnCompLengths(n)= length(CC.PixelIdxList{n});
end
clear n;

[max1, indx1] = max(ConnCompLengths); %find longest components and its index
ConnCompLengths(indx1) = 0; % delete that max, to find the second maximum
[max2, indx2] = max(ConnCompLengths); %find second largest components

[r1,tmp] = ind2sub(size(I),CC.PixelIdxList{indx2}); %find the row/col index
r2 = max(r1); % top row for the first component
r1 = min(r1); % find the minimum row in one component of sobel

newI((1:r1),:) = 0; %Top rows will be deleted, to limit the area of interest between r1 and r2


[rowsI,columnsI] = size(I);

boundaries_hotpixel = zeros(size(I)); %empty image

%limit search area between r1 and r2, and find the boundaries in newI using
%the two loops below. 
for n = r1:r2
    if (sum(newI(n,:))~=0) %the sum per row
        loc = find(newI(n,:)>0); %find non-zero values in that row
        boundaries_hotpixel(n,loc(1)) = 2^16; %take first non-zero pixel in each row, and set it to 1
        boundaries_hotpixel(n,loc(end)) = 2^16; %take the last non-zero pixel in each row
    end   
end

for n = 1:columnsI
    if (sum(newI((r1:r2),n))~=0) %the sum per col
        loc = find(newI((r1:r2),n)>0); %find non-zero values in that col
        boundaries_hotpixel((loc(1)+r1-1),n) = 2^16; %take first non-zero pixel in each col, and set it to 1
    end   
end
clear n;
% figure, imshow(boundaries), title('Boundaries with HotPixelFinder')

[row1, col1] = find(sum(boundaries_hotpixel>0));

for f=1:length(col1) %go column by column 
    if sum(boundaries_canny(:,f))>0 %if there are pixels in that column in canny
       boundaries_hotpixel(:,f) = 0; %delete that column in boundary_hotpixel 
    end    
end
clear f


loca = find(boundaries_hotpixel>0); %find locations of all non-zero in boundaries
boundaries_hotpxcanny(loca) = 2^16; %Keep Hotpixel edges as well

figure, imshow(boundaries_hotpxcanny), title('Canny edges AND HotPixelFinder')


figure, imshow(I,[]), title('Red: HotPixel finder AND Canny combined')
% red on top on figure
red = cat(3, ones(size(I)), zeros(size(I)), zeros(size(I))); %red has RGB value 1 0 0
hold on 
h = imshow(red); 
hold off 
% Use our diff1 as the AlphaData for the solid red image. 
set(h, 'AlphaData', boundaries_hotpxcanny) 

%---------------------------------------------------------------------


toc

kefpart_2 %move on to part 2 of code
