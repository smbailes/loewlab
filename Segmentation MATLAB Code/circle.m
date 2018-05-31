% Circle Test Code for Hough Transform

% This code aims to see if ellipse code will give the same parameters for 
% major and minor axes if it is asked to find an ellipse that is a perfect
% circle

% Last Updated by Zainab Mahmood 11/9/17 at 12:14PM


%% Part 1: Create circle

clear all;
clc;
close all;


% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = 640;
imageSizeY = 480;
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = 320;
centerY = 240;
radius = 100;
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;
% circlePixels is a 2D logical array.
% Now, display it.
image(circlePixels) ;
colormap([0 0 0; 1 1 1]);
title('Binary image of a circle');

%% Part 2: Segment Circle Using Hough Transform

edgecanny = edge(circlePixels,'canny');    %find canny edges
figure,imshow(edgecanny);
title('Canny Edges');
    
figure;
imshow(circlePixels);
    
% set parameters for ellipse detection
paramsr.minMajorAxis = 210;
paramsr.maxMajorAxis = 500;
paramsr.numBest = 12; %draws 12 ellipses
paramsr.rotation = 45; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
paramsr.rotationSpan = 35;
paramsr.randomize = 7; %randomization component that may reduce changing of
%ellipses

% note that the edge (or gradient) image is used
bestFitsr = ellipseDetection(edgecanny, paramsr);
fprintf('Output %d best fits.\n', size(bestFitsr,1));

%ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 

%takes the information that was found of the ellipses and draws them;also
%keeping the information for each ellipse in a cell in qr(and later ql for those):
qr{1,length(bestFitsr)}=0;
for n=1:length(bestFitsr)
    qr{n} = ellipse(bestFitsr(n,3),bestFitsr(n,4),bestFitsr(n,5)*pi/180,bestFitsr(n,1),bestFitsr(n,2),'k');
end


%ELLIPSES TO PIXELS
[img_y, img_x] = size(circlePixels);
ellipses=zeros(img_y,img_x);

%draw in circle pixel by pixel
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
        ellipses(ye2(a,d),xe2(a,d)) = 1;      %fill in 1's wherever there is a point in the ellipse
    end
end

% show red on top of figure
red = cat(3, ones(size(circlePixels)), zeros(size(circlePixels)), zeros(size(circlePixels))); %red has RGB value 1 0 0
hold on 
displ = imshow(red); 
hold off 
set(displ, 'AlphaData', ellipses)
 

