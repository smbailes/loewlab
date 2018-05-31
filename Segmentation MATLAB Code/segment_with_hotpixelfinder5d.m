% Segment the image with HotPixelFinder2_uint16
% Last modified: 05/01/2017 11:42 AM
% *****LOOK AT DIFFERENCES IN PART 5******
% This copy of segment_with_hotpixelfinder.m does the following:
%       1) reflect the breast boundary if necessary
%       2) It connects the middle regions and to the side boundaries
%           to make a one-piece boundary 
%       3) is will connect the middle boundaries BEFORE adding SOBEL to the
%       image
%       4) This finds if canny edge and hotpixel are in the same column,
%       delete hotpixel and keeps canny
% It will remove rows above the minimum (top) row in the sobel boundaries,
% hence, deleting any extra lines that would affect the results and limit
% the edges to only what's below that threshold row. It'll be removed from
% the hotpixelfinder image (newI)


%% Part 1: Hot pixel finder
% Ask the user for a percentage, and find the top p% of pixels 
% Create a new binary image of that top p% region.

clear all;
clear;
close all;
clc;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
I = imread(ptID); %open the image, keeping it in 16-bits

bins = 2^16; %insert image bits here
[N,binlocation] = imhist(I,bins); %each count will has its own bin

perc = input('What is your desired percentage? '); %insert user desired percentage
hotPix = 0; 
lowerBound = .009*perc*numel(I); % numel(I): the number of pixels in the image
upperBound = .011*perc*numel(I);
global imTitle;
imTitle = sprintf('Highest %d percent of %s',perc, ptID);
total = 0; 
for j = 1:numel(N);  %Shouldn't it be 0 to numel(N)
    total = total+ N(bins-j); %Tells us where the 1 percent of pixels start from but some of the bins still have zero in them
    if (total >= lowerBound)  
       break
    end 
end


k = 0;
for i = bins-j:bins; %Keep all values that say 256 as 256. Matlab does not have a bin called 0 so there is an extra bin.
    if (N(i)~=0); %If the bin does not equal zero then move on to the next and add it to the hotPix list 
     k = k+1; 
     hotPix(k) = i;  %This contains the 1 percent of pixels that do not contain zero
     
    end 
end

newI = zeros(size(I)); %new image and make the entire image black
for m = 1:k;
    o = find(I==hotPix(m)); %This looks for the hotPix in the image
    newI(o) = bins; %Make any value of the bightest 1 percent white
end
figure, imshow(newI)
title(imTitle)

clear lowerBound upperBound hotPix bins binlocation total;
%% Part 2.5: add in ellipse code


e1 = edge(I,'canny');%finds canny edge
CCANNY = bwconncomp(e1);
figure, imshow(I,[])
title('Image 1: Original Image')
figure, imshow(e1) %displays canny edges of p
title('Image 2: Canny Image, before small lines removed')
pxlength = 10; %pixels to be removed
boundededges = bwareaopen(e1, pxlength); %Remove short lines
figure, imshow(boundededges)
title('Image 3: Canny Image, small lines removed')
smaller = bwmorph(boundededges,'skel',Inf); %move the canny edges down to skeleton
figure, imshow(smaller)
title('Image 4: Skeleton Canny Image, to find ellipses from')

figure;
imshow(I,[]); %produce an image to overlay the ellipses onto
title('Final Image with Ellipses')

%left side
% override some default parameters
paramsl.minMajorAxis = 200;
paramsl.maxMajorAxis = 300;
paramsl.numBest = 6; %draws 6 ellipses
paramsl.rotation = 135; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
paramsl.rotationSpan = 35;

% note that the edge (or gradient) image is used
bestFitsl = ellipseDetection(smaller, paramsl);
fprintf('Output %d best fits.\n', size(bestFitsl,1));

%ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
for n=1:length(bestFitsl)
    ql{n} = ellipse(bestFitsl(n,3),bestFitsl(n,4),bestFitsl(n,5)*pi/180,bestFitsl(n,1),bestFitsl(n,2),'k');
end
saveellipsesl = ellipse(bestFitsl(:,3),bestFitsl(:,4),bestFitsl(:,5)*pi/180,bestFitsl(:,1),bestFitsl(:,2),'k');
%the above is trying to use a handle to keep the ellipses....not much from
%this yet.
%!! Will be changing and separating the two sides of the image for the two
%breasts !!

%right side
paramsr.minMajorAxis = 300;
paramsr.maxMajorAxis = 400;
paramsr.numBest = 6; %draws 6 ellipses
paramsr.rotation = 45; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
paramsr.rotationSpan = 35;

% note that the edge (or gradient) image is used
bestFitsr = ellipseDetection(smaller, paramsr);
fprintf('Output %d best fits.\n', size(bestFitsr,1));

%ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 

for n=1:length(bestFitsr)
    qr{n} = ellipse(bestFitsr(n,3),bestFitsr(n,4),bestFitsr(n,5)*pi/180,bestFitsr(n,1),bestFitsr(n,2),'k');
end
saveellipsesr = ellipse(bestFitsr(:,3),bestFitsr(:,4),bestFitsr(:,5)*pi/180,bestFitsr(:,1),bestFitsr(:,2),'k');



%ellipse(semimajor axis of ra,semimajor axis of radius rb, semimajor axis
%of ang, center point x0, center pointy0, 'color')
%
% **The length of ra, rb, and ang should be the same.**

%the following is for finding the Hough circles, and then overlaying them
rmin = 75;
rmax = 125;
[centers, radii] = imfindcircles(smaller,[rmin rmax],'Sensitivity',0.97); %A higher 'Sensitivity' value sets the detection threshold lower and leads to detecting more circles.
h = viscircles(centers,radii,'LineWidth',1);


%% Part 2: find the boundaries of hot regions in newI.
% Locate the one-pixel thick boundary of the hot regions in newI:
%       - The result can be used for Part 3, 
%       or the entire hot region can be used instead of just its boundaries.
%       - I chose to use the entire region, because it gave better results
%       and you can maximize the overlap area (See Part 3)
% *****  This part MUST be improved *****

% ** Suggested improvement: limit the search to the area between max and
% min row of sobel boundaries:
%       - sobel boundaries detect the two body side boundaries perfectly.

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

newI((1:r1),:) = 0; %Top rows will be deleted, to limit the area of interest between r1 and r2


[r c] = size(I);

boundaries_hotpixel = zeros(size(I)); %empty image

%limit search area between r1 and r2, and find the boundaries in newI using
%the two loops below. 
for n = r1:r2
    if (sum(newI(n,:))~=0) %the sum per row
        loc = find(newI(n,:)>0); %find non-zero values in that row
        boundaries_hotpixel(n,loc(1)) = 2^16; %take first non-zero pixel in each row, and set it to 1
        boundaries_hotpixel(n,loc(end)) = 2^16; %take the last non-zero pixel in each row
    end   
    n = n + 1;
end

for n = 1:c
    if (sum(newI((r1:r2),n))~=0) %the sum per col
        loc = find(newI((r1:r2),n)>0); %find non-zero values in that col
        boundaries_hotpixel((loc(1)+r1-1),n) = 2^16; %take first non-zero pixel in each col, and set it to 1
    end   
    n = n + 1;
end
clear n;
% figure, imshow(boundaries), title('Boundaries with HotPixelFinder')


%% Part 3: Match with edge detection techniques, and only keep the canny boundaries
% Slowest part in the code
% Find the connected components in I_canny, then look if any of these
% components overlap with the hotregions located in newI. 
% If at least one pixel in a canny component fall overlaps with the
% hotregion in newI, then store the ENTIRE canny component

I_canny = edge(I,'canny');

% Connecting the canny boundaries could be helpful eventually, but for now
% it didn't work --- I suggest playing around with this.
% se = strel('disk',1); %Create a Morphological structuring element, you change the shape used and diameter
% I_canny= imclose(I_canny,se); % Morphologically close image using the se element
% figure,imshowpair(I_canny,I_sobel,'montage'),title('canny and sobel')

CC = bwconncomp(I_canny); %find connected components in canny


boundaries_canny = zeros(size(I));
loc = find(newI>0); %match with the entire area in the p% from above part
for x = 1:length(loc) %go index by index
    for y = 1:CC.NumObjects %go through each connected component in CC
        if find(CC.PixelIdxList{y}==loc(x)) %check if location x is found in any of the connected components
            %curly brackets access the contents of the the cell
            boundaries_canny(CC.PixelIdxList{y}) = 2^16; %if yes, then keep the ENTIRE component and store it in boundaries2
        end
    end
end
clear x y loc;
% figure, imshow(boundaries2), title('Canny edges using HotPixelFinder as a reference')

%% Part 4: Display the detected edges on the original image
% For the future I suggest writing this code in a function file, and just
% call the function in here to increase efficiency. Because this part was 
% repeated at least 3 times.

imgOG = imread_irfanview(ptID); %for display purposes only
figure, imshow(imgOG), title('Red: HotPixel finder. Blue: Canny after using HotPixel as a reference')

% red on top on imgOG
red = cat(3, ones(size(imgOG)), zeros(size(imgOG)), zeros(size(imgOG))); %red has RGB value 1 0 0
hold on 
h = imshow(red); 
hold off 
% Use our diff1 as the AlphaData for the solid red image. 
set(h, 'AlphaData', boundaries_hotpixel) 

% Overlay the detected region on the original image, and compare with
% HotPixel Edges

blue = cat(3, zeros(size(imgOG)), zeros(size(imgOG)), ones(size(imgOG))); %blue has RGB value of 0 0 1
hold on 
h2 = imshow(blue); 
hold off 
set(h2, 'AlphaData', boundaries_canny) 

%% Part 5: Combine the edges from HotPixel finder AND Canny

boundaries_hotpxcanny = boundaries_canny;

[row1, col1] = find(sum(boundaries_hotpixel>0));

for n=1:length(col1) %go column by column 
    if sum(boundaries_canny(:,n))>0 %if there are pixels in that column in canny
       boundaries_hotpixel(:,n) = 0; %delete that column in boundary_hotpixel 
    end    
end

loc = find(boundaries_hotpixel>0); %find locations of all non-zero in boundaries
boundaries_hotpxcanny(loc) = 2^16; %Keep Hotpixel edges as well
clear x y loc;

%figure, imshow(boundaries3), title('Canny edges AND HotPixelFinder')

% imgOG = imread_irfanview(ptID);
% figure, imshow(imgOG), title('Red: HotPixel finder AND Canny combined')
% % red on top on imgOG
% red = cat(3, ones(size(imgOG)), zeros(size(imgOG)), zeros(size(imgOG))); %red has RGB value 1 0 0
% hold on 
% h = imshow(red); 
% hold off 
% % Use our diff1 as the AlphaData for the solid red image. 
% set(h, 'AlphaData', boundaries3) 

%% Part 6: Refine the combined edges resulting from Part 5
% Merge the edges from hotpixel and canny.

% Dilation is an option, haven't implemented it, but it's worth playing
% with
%    - an isotropic dilation step: are "almost connected" if they are within 4 
%    pixel units of distance from each other. 
% boundaries3 = bwdist(boundaries3) <= 2;
% figure, imshow(boundaries3)

se = strel('disk',4); %Create a Morphological structuring element, you change the shape used and diameter
boundaries_hotpxcanny = imclose(boundaries_hotpxcanny,se); % Morphologically close image using the se element

%% Part 7: Add the Sobel Part temporarly (The two body sides)
% This is to get the top and min rows in order to limit the reflection 
% and search region for part 8. This way the top part of the image above
% the armpits will be cropped (hopefully)

I_sobel = edge(I,'sobel');
se = strel('disk',2); %Create a Morphological structuring element, you change the shape used and diameter
I_sobel= imclose(I_sobel,se); % Morphologically close image using the se element

CC = bwconncomp(I_sobel);

%Find the length of each conncomp in CC and store them in an array
for n = 1:CC.NumObjects
    ConnCompLengths(n)= length(CC.PixelIdxList{n});
end
clear n;

[max1, indx1] = max(ConnCompLengths); %find longest components and its index
ConnCompLengths(indx1) = 0; % delete that max, to find the second maximum
[max2, indx2] = max(ConnCompLengths); %find second largest components

[r1,tmp] = ind2sub(size(I),CC.PixelIdxList{indx2});
r2 = max(r1); % top row for the first component
r1 = min(r1); % find the bottom row

% Remove any connected components (CANNY+HOTPIX) that are above the top row
% of sobel ---- *** This does not work for all cases, especially if some of
% the canny boundaries of interest extend beyond the armpits (Look at
% patient IRST002)
%
% boundaries3((1:r1),:) = 0; %Top rows will be deleted

%Use indx1 and indx2 in the next part.

%% Part 8: Reflect one of the boundaries to the other side if necessary
%NOTE: boundaries_hotpxcanny DOESN'T have sobel edges yet

figure, imshow(boundaries_hotpxcanny);

%Remove short lines if necessary, you don't want extra lines to be
%reflected or cause errors. The more lines you get rid of the better.
while(1) %Keep asking the user to remove short segments, until they break 
    inpt = input('Do you want to remove short segments [Y/N]? ','s');
    if inpt=='y'
        pxlength = input('What''s the max conncomp length you want to remove? '); %Enter number of pixels
        boundaries_hotpxcanny = bwareaopen(boundaries_hotpxcanny, pxlength); %Remove short lines
        imshow(boundaries_hotpxcanny)
    else
        break;
    end
end
clear inpt pxlenth;

%Specify region to look at.. between the two sobel boundaries
% Must find the two columns between the sobel boundaries that would limit
% the search for the middle region.
%   - So far, we have the indeces of the two sobel components
%   - If indx1 is on the left side of the body and indx2 on the right, 
%   decide which one comes first and determine the two columns we need (c1
%   and c2).
if CC.PixelIdxList{indx1}(end) < CC.PixelIdxList{indx2}(1)
    [tmp,c1] = ind2sub(size(I),CC.PixelIdxList{indx1}(end))
    [tmp,c2] = ind2sub(size(I),CC.PixelIdxList{indx2}(1))
end

% If indx1 is on the right side of the body and indx2 on the left, 
% decide which one comes first
if CC.PixelIdxList{indx1}(1) > CC.PixelIdxList{indx2}(end)
    [tmp,c1] = ind2sub(size(I),CC.PixelIdxList{indx2}(end))
    [tmp,c2] = ind2sub(size(I),CC.PixelIdxList{indx1}(1))
end
clear tmp;

% find middle point between c1 and c2 --- This will be the reflection
% mirror
% ** I suggest finding a better way to locate the middle column
middleCol = round((c1+c2)/2);

%Method A: Summing columns
sumCol = sum(boundaries_hotpxcanny((r1:r2),(c1:c2))); %sum columns between c1 and c2 only
% Keep between r1 and r2 to remove excess pixels from the top of the image.
%figure, plot((c1:c2),sumCol), title('Sum of Columns between the 2 sobel boundaries')

zeroIndx = find(sumCol==0); %zeroIndx has value of a column num, 
% but on the original image it should be that value PLUS c1
for n = 1:length(zeroIndx)
    if (zeroIndx(n)+c1) > middleCol   %The zero index is on the right breast(or patient's left)
       if (sumCol(zeroIndx(n) + c1 - 2*(zeroIndx(n) + c1 - middleCol)) ~=0) %Find mirror index on the left side
           %if the the sum of the corrosponding column on the other side is 
           % NOT zero, meaning there is a a boundary value somewhere in 
           % that column
           
           %Copy the column on the left to the column on the right
           boundaries_hotpxcanny(:,(zeroIndx(n)+c1)) = boundaries_hotpxcanny(:,(zeroIndx(n) + c1 -  2*(zeroIndx(n) + c1 - middleCol)));
       end
    end
    if (zeroIndx(n)+c1) < middleCol  % zero index is on the left breast
         if (sumCol(zeroIndx(n) + 2*(middleCol - zeroIndx(n) - c1)) ~=0) %Find mirror index on the right side
           %if the the sum of the mirror column is NOT zero, meaning there
           %is a a boundary value somewhere in that column
           
           %Copy the column on the right to the column on the left
           boundaries_hotpxcanny(:,(zeroIndx(n)+c1)) = boundaries_hotpxcanny(:,(zeroIndx(n) +  2*(middleCol - zeroIndx(n)) - c1));
         end
    end
end

% Thought: Is there a way to check with canny?

%figure, imshow(boundaries3)
se = strel('disk',2); %Create a Morphological structuring element, you change the shape used and diameter
boundaries_hotpxcanny = imclose(boundaries_hotpxcanny,se); % Morphologically close image using the se element

%overlay on the original image
imgOG = imread_irfanview(ptID);
figure, imshow(imgOG), title('Reflecting one boundary over the other to fill in spaces')
% red on top on imgOG
red = cat(3, ones(size(imgOG)), zeros(size(imgOG)), zeros(size(imgOG))); %red has RGB value 1 0 0
hold on 
h = imshow(red); 
hold off 
% Use our diff1 as the AlphaData for the solid red image. 
set(h, 'AlphaData',boundaries_hotpxcanny) 

clear i j k m n N o r c max1 max2 imTitle perc ;

%% Part 9: Connect the breast boundaries in the middle area only
%   Still doesn't have sobel

%Reduce the boundaries to one-pixel thick lines, without disconnecting them
boundaries4 = bwmorph(boundaries_hotpxcanny,'skel',Inf);

figure, imshow(boundaries4), title('Morphological operation: Skeleton')
while(1) %Keep asking the user to remove short segments, until they break out of loop
    inpt = input('Do you want to remove short segments [Y/N]? ','s');
    if inpt=='y'
        pxlength = input('What''s the max conncomp length you want to remove? ');
        boundaries4 = bwareaopen(boundaries4, pxlength); %Remove short lines
        imshow(boundaries4)
    else
        break;
    end
end
clear inpt pxlenth;

CC2 = bwconncomp(boundaries4);
newboundaries4 = boundaries4;

%Idea is to connect the last pixel from conncomp 1 with the first pixel on
%conncomp 2. The first and last pixels are determined by the coloums. Last
%pixel is the pixel with the max column index, and first pixel is with the
%min column index on the connected component. Remember that the columns
%indeces are numbered from left to right of the image.

for n = 1:CC2.NumObjects - 1 %the number of lines in the middle region of the patient
    %Store all row and col values of component n and the component after in
    % x1,y1, x2, y2
    [x1, y1] = ind2sub(size(newboundaries4),CC2.PixelIdxList{n});
    [x2, y2] = ind2sub(size(newboundaries4),CC2.PixelIdxList{n+1});
    
    [yy1, ind] = max(y1); %find the max col in component n
    xx1 = x1(ind); % The corrosponding row value for max col
    
    [yy2, ind] = min(y2); %find the min col in component n+1
    xx2 = x2(ind); % The corrosponding row value for min col
    
    %Draw a line between the two points (xx1,yy1) and (xx2,yy2) and insert
    %it in newboundaries4
    shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White');
    newboundaries4 = step(shapeInserter, newboundaries4, uint16([yy1 xx1 yy2 xx2]));
    %figure, imshow(newboundaries4), title('After step shapeinserter');
    
end
clear xx2 xx1 yy1 yy2 y1 y2 x1 x2;

figure, imshow(newboundaries4), title('Connecting the middle region')


%% Part 10: Add the Sobel Part - Include the two body sides.

I_sobel = edge(I,'sobel');
se = strel('disk',2); %Create a Morphological structuring element, you change the shape used and diameter
I_sobel= imclose(I_sobel,se); % Morphologically close image using the se element

CC = bwconncomp(I_sobel);

%Find the length of each conncomp in CC and store them in an array
for n = 1:CC.NumObjects
    ConnCompLengths(n)= length(CC.PixelIdxList{n});
end
clear n;

[max1, indx1] = max(ConnCompLengths); %find longest component and its index
ConnCompLengths(indx1) = 0; % delete that max, to find the second maximum
[max2, indx2] = max(ConnCompLengths); %find second largest component

%newboundaries4b = newboundaries4; %To keep original intact, and doesn't have sobel

%threshold the 2 longest edges/components from sobel to the image
newboundaries4(CC.PixelIdxList{indx1}) = 2^16 ; %This has the SOBEL and CANNY
newboundaries4(CC.PixelIdxList{indx2}) = 2^16;

%figure, imshow(boundaries3), title('Boundaries with Sobel combined')
%boundaries4 = bwareaopen(boundaries3,100); %Remove short lines
%boundaries3b = bwareaopen(boundaries3,100); 

%overlay on the original image
imgOG = imread_irfanview(ptID);
figure, imshow(imgOG), title('HotPixel finder AND Canny AND Sobel combined after refining')
% red on top on imgOG
red = cat(3, ones(size(imgOG)), zeros(size(imgOG)), zeros(size(imgOG))); %red has RGB value 1 0 0
hold on 
h = imshow(red); 
hold off 
% Use our diff1 as the AlphaData for the solid red image. 
set(h, 'AlphaData',newboundaries4) 

% clear CC max1 max2 indx1 indx2 n;
%% Part 11: Connect the breast boundaries to the sides

%Reduce the boundaries to one-pixel thick lines, without disconnecting them
newboundaries4 = bwmorph(newboundaries4,'skel',Inf);
figure, imshow(newboundaries4), title('Morphological operation: Skeleton')

while(1) %Keep asking the user to remove short segments, until they break 
    inpt = input('Do you want to remove short segments [Y/N]? ','s');
    if inpt=='y'
        pxlength = input('What''s the max conncomp length you want to remove? ');
        newboundaries4 = bwareaopen(newboundaries4, pxlength); %Remove short lines
        imshow(newboundaries4)
    else
        break;
    end
end
clear inpt pxlenth;

% Now connect the middle boundary to the two body sides:
CC2 = bwconncomp(newboundaries4);

% ***** CONNECT THE LEFT SIDE TO THE MIDDLE BOUNDARY:
%Store all row and col values of component n and the component after in
% x1,y1, x2, y2
[x1, y1] = ind2sub(size(newboundaries4),CC2.PixelIdxList{1});
[x2, y2] = ind2sub(size(newboundaries4),CC2.PixelIdxList{2});

[yy2, ind] = min(y2); %find the min col in component 2
xx2 = x2(ind); % The corrosponding row value for min col
% ^^ this is the first point

%Calculate the distance between xx2,yy2 and each point on the sobel side
%edge, then find the point with the shortest distance.

% **** LEFT BREAST (or patient's right)
distance = [];
for n=1:length(x1)
   xx1 = x1(n); %Take coordinates of one pixel at a time from the sobel boundary
   yy1 = y1(n);
    %Find distance between the two points and store it in distance
   distance = [distance pdist([xx1,yy1;xx2,yy2])]; %Euclidean distance
end
Distances = [x1 y1 distance']; %Arrange each point with its distance in a matrix
[minDist loc] = min(Distances(:,3)) %find the smallest distance and its location in the matrix

% Find the second point
xx1 = Distances(loc,1); %from col 1 at loc
yy1 = Distances(loc,2); %from col 2 at loc
   
%Draw a line between the two points (xx1,yy1) and (xx2,yy2)
shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White');
newboundaries4 = step(shapeInserter, newboundaries4, uint16([yy1 xx1 yy2 xx2]));

% ***** CONNECT THE Right SIDE TO THE MIDDLE BOUNDARY:
%Store all row and col values of component n and the component after in
% x1,y1, x2, y2
[x1, y1] = ind2sub(size(newboundaries4),CC2.PixelIdxList{end});
[x2, y2] = ind2sub(size(newboundaries4),CC2.PixelIdxList{end-1});

[yy2, ind] = max(y2); %find the max col in component 2
xx2 = x2(ind); % The corrosponding row value for min col
% ^^ this is the first point

% **** RIGHT BREAST (or patient's left)
distance = [];
for n=1:length(x1)
   xx1 = x1(n);
   yy1 = y1(n);
   distance = [distance pdist([xx1,yy1;xx2,yy2])];
end
Distances = [x1 y1 distance'];
[minDist loc] = min(Distances(:,3))

% Find the second point
xx1 = Distances(loc,1); %from col 1 at loc
yy1 = Distances(loc,2); %from col 2 at loc
   
%Draw a line between the two points (xx1,yy1) and (xx2,yy2)
shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White');
newboundaries4 = step(shapeInserter, newboundaries4, uint16([yy1 xx1 yy2 xx2]));
figure, imshow(newboundaries4), title('Final image with all boundaries connected')
   

%overlay on the original image
imgOG = imread_irfanview(ptID);
figure, imshow(imgOG), title('All edges connected (hopefully!)')
% red on top on imgOG
red = cat(3, ones(size(imgOG)), zeros(size(imgOG)), zeros(size(imgOG))); %red has RGB value 1 0 0
hold on 
h = imshow(red); 
hold off 
% Use our diff1 as the AlphaData for the solid red image. 
set(h, 'AlphaData',newboundaries4) 

disp('End of Program')