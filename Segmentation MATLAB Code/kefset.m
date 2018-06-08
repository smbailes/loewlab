
close all;
clear;
clc;

%% Part 1 Intial Loading of Images

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
dir = uigetdir; 
I = imread([dir '\' ptID]); 

figure, imshow(I,[]) %to help decide if it should be cropped or not
title('Original Image')

% cropopt=input('Does a crop need to occur?? y/n: ','s');
% if cropopt == 'y'
%     I=imcrop(I,[]); %cropping, if necessary
% end

I=getMatrixOutliers(I);
figure, imshow(I,[]) %to help decide if it should be cropped or not
title('Outliers Removed')

close Figure 1
%perc = input('What is your desired percentage? '); 
perc = 5; %top 5% of pixels used

edgecanny = edge(I,'canny');
edgecanny=bwareaopen(edgecanny,10); %removes very small edge lines

figure,imshow(edgecanny)
title('Canny edges');

%% Part 2, Ellipse Detection

figure;
imshow(I,[]); %produce an image to overlay the ellipses onto
title('Image with Ellipses')

in = input('Is the breast small or large? Enter s/l: ','s');
% if in == 's'
%         %RIGHT SIDE
%     % override some default parameters
%     paramsr.minMajorAxis = 150;
%     paramsr.maxMajorAxis = 300;
%     paramsr.numBest = 12; %draws 12 ellipses
%     paramsr.rotation = 45; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
%     paramsr.rotationSpan = 35;
%     paramsr.randomize = 7; %randomization component that may reduce changing of
%     %ellipses
% 
%     % note that the edge (or gradient) image is used
%     bestFitsr = ellipseDetection(edgecanny, paramsr);
%     fprintf('Output %d best fits.\n', size(bestFitsr,1));
% 
%     %ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
% 
%     %takes the information that was found of the ellipses and draws them;also
%     %keeping the information for each ellipse in a cell in qr(and later ql for those):
%     qr{1,length(bestFitsr)}=0;
%     for n=1:length(bestFitsr)
%         qr{n} = ellipse(bestFitsr(n,3),bestFitsr(n,4),bestFitsr(n,5)*pi/180,bestFitsr(n,1),bestFitsr(n,2),'k');
%     end
%     %overriding parameters:
%     paramsl.minMajorAxis = 150;
%     paramsl.maxMajorAxis = 300;
%     paramsl.numBest = 12; %draws 12 ellipses
%     paramsl.rotation = 135; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
%     paramsl.rotationSpan = 35;
%     paramsl.randomize = 7; %randomization component that may reduce changing of
%     %ellipses
%     
%     
%     %LEFT SIDE
%     bestFitsl = ellipseDetection(edgecanny, paramsl);
%     fprintf('Output %d best fits.\n', size(bestFitsl,1));
%     
%     %ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
%     ql{1,length(bestFitsl)}=0;
%     for n=1:length(bestFitsl)
%         ql{n} = ellipse(bestFitsl(n,3),bestFitsl(n,4),bestFitsl(n,5)*pi/180,bestFitsl(n,1),bestFitsl(n,2),'k');
%     end
%     
if in == 'l'
    %RIGHT SIDE
    % override some default parameters
    paramsr.minMajorAxis = 350;
    paramsr.maxMajorAxis = 700;
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
    %overriding parameters:
    paramsl.minMajorAxis = 350;
    paramsl.maxMajorAxis = 700;
    paramsl.numBest = 12; %draws 12 ellipses
    paramsl.rotation = 135; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
    paramsl.rotationSpan = 35;
    paramsl.randomize = 7; %randomization component that may reduce changing of
    %ellipses


    %LEFT SIDE
    bestFitsl = ellipseDetection(edgecanny, paramsl);
    fprintf('Output %d best fits.\n', size(bestFitsl,1));

    %ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
    ql{1,length(bestFitsl)}=0;
    for n=1:length(bestFitsl)
        ql{n} = ellipse(bestFitsl(n,3),bestFitsl(n,4),bestFitsl(n,5)*pi/180,bestFitsl(n,1),bestFitsl(n,2),'k');
    end


%ELLIPSES TO PIXELS
[img_y, img_x] = size(I);
ellipses=zeros(img_y,img_x);

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
        ellipses(ye1(a,d),xe1(a,d)) = 1;      %fill in 1's wherever there is a point in the ellipse
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
        ellipses(ye2(a,d),xe2(a,d)) = 1;      %fill in 1's wherever there is a point in the ellipse
    end
end

[checky,checkx]=size(ellipses);
if checkx > img_x
    ellipses=ellipses(1:img_y,1:img_x);
end
if checky > img_y
    ellipses=ellipses(1:img_y,1:img_x);
end

sf = strel('disk',2); %Create a Morphological structuring element, you change the shape used and diameter
ellipses= imclose(ellipses,sf); 

ellipses=bwmorph(ellipses,'clean');

ellipses(1:round(img_y/4),:)=0;

if in=='s'
    ellipses(:,1:round(img_x/3))=0;
    ellipses(:,round(2*(img_x/3)):end)=0;
end


figure, imshow(I,[]), title('Ellipses')
% red on top on figure
red = cat(3, ones(size(I)), zeros(size(I)), zeros(size(I))); %red has RGB value 1 0 0
hold on 
displ = imshow(red); 
hold off 
set(displ, 'AlphaData', ellipses)

end

%% Part 3, Hot Pixel Finder

bins = 2^16; %insert image bits here
[N,binlocation] = imhist(I,bins); %each count will has its own bin


lowerBound = .009*perc*numel(I); % numel(I): the number of pixels in the image
upperBound = .011*perc*numel(I);

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
title('Hot Pixel');

%% Part 4: Small/Large 
figure, imshow(I,[]), title('Bound Detection')
%hold on; 
fprintf('Select Upper Bound \n');
[Xup,Yup] = ginput(1);
%plot(Yup, :)
%hold on;
fprintf('Select Lower Bound \n');
[Xlo,Ylo] = ginput(1);
%plot(Ylo,:,'g')
%hold on; 
fprintf('Select Left Bound \n');
[Xleft,Yleft] = ginput(1);
%plot(:,Xleft,'g')
%hold on; 
fprintf('Select Right Bound \n');
[Xright,Yright] = ginput(1);
%plot(:,Xright,'g')
%hold off;
if in == 's'
     kefsmall1
elseif in == 'l'
     keflarge1
end


