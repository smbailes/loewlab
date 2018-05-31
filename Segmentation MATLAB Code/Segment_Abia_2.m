% Segmentation Algorithm 2 
% Author : Abia Khan
% Edited by: Nada Kamona
% Last revised by Author: 12/13/2016
% Last Edited by editor: 2/4/2017  7:00 AM

% The segmentation method is below. 

% This is an attempt to make the code works for IRST001, IST003, and
% IRST005.
% *** NOT working yet ***

% Note: Test with and without removing background, see how that affects
% edge detection with canny.

%% Step 1: Preprocessing: Input Image, Clear background, Revised egde detection algorithm
clear;
close all;
clc;
% Read Image
sa=imread_irfanview('1799 - P5.tif');

% Crop Image
sa=imcrop(sa);

%Clear Background
%thresholdvalue = graythresh(sa)*2^16; %Clearing the background : Otsu's method

%Nada: The method above is removing pixels of interest in the breast
%region, and it's making the canny edges less apparent. 

%Nada: New threshold method
% thresholdvalue = 7000; %Clearing the background

% [row_left,col_left]=size(sa);
% 
% for k = 1:row_left  %Nada: Can replace this loop with find function
% 
%     for u = 1:col_left
% 
%         if sa(k,u)<thresholdvalue
% 
%             sa(k,u)=0; 
% 
%         end
% 
%     end
% 
% end
% imshow(sa)

edge=edge(sa,'canny'); % Edge Detection
imshow(edge)
title('Canny only')

se = strel('disk',5); %Create a Morphological structuring element, you change the shape used and diameter
bwImgs = imclose(edge,se); % Morphologically close image using the se element
figure, imshow(bwImgs)
title('closing lines first time')

bwImgs = bwareaopen(bwImgs,200); % removes 200 pixes with a connectivity of 8 in the image. This gets rid of almost all of the vessels/spurious lines in the breast.
figure
imshow(bwImgs)
title('Removing short lines')

%bwImgs = bwmorph(bwImgs,'close',inf); %morphological closing (dilation followed by erosion)

se2 = strel('disk',5); %Create a Morphological structuring element, you change the shape used and diameter
bwImgs = imclose(bwImgs,se2); % Morphologically close image using the se element

figure,imshow(bwImgs);
title('Filling in between the lines using morphological element')

% bwImgs = bwareaopen(bwImgs,100); % removes 200 pixes with a connectivity of 8 in the image. This gets rid of almost all of the vessels/spurious lines in the breast.
% figure
% imshow(bwImgs)
% title('Removing short lines round 2')

bounds = bwboundaries(bwImgs); % Computes the boundaries in the image


%% Step 2: Display right and left breast boundaries

boundary = bounds{2}; %gives the properties of the second boundary, which is the right breast
%  figure
%  plot(boundary(:,2), boundary(:,1),'r'); 

row_right=boundary(:,1); % rows of right breast boundary

col_right=boundary(:,2); % columns of right breast boundary

boundaryl = bounds{3};  %Nada: Should be the left breast boundry

row_left=boundaryl(:,1); % rows of left breast boundary

col_left=boundaryl(:,2); % columns of left breast boundary


 %% Perform Segmentation
%  
% figure
% imshow(sa)
y1=get(gca,'ylim'); % the y limit of the image

w_a= min(col_right); %Lowest y coordinate for right boundary
w_b=max(col_right); %Highest y coordinate for right boundary

%qw,qe, qr, qt are the coordinates of the top two pixels on the left and
%right breast boundaries.
qw = max(col_right); %qw, qe are the last pixel on the left breast (highest point in the middle region)
qe = min(row_right);

qr = min(row_left); %qr, qt are the first pixel on the left breast (highest point in the middle region)
qt = min(col_left);


y_i= [qe qr]; % Create coordinates
x_i= [qw qt];  % Create coordinates


gh = [qw:qt]'; %Create rows between the two breast boundaries (middle chest region)
%gh=gh';
fg=[qr:qe]';
%fg=fg';
 

f = polyfit(x_i,y_i,1);%Draw best fit line between the region between the boundaries
f=(f(1));

bv = (qr-(f*qt)); % Perform the slope
vb= (f*qw)+bv;

bn = (f*gh)+bv; %Line of best fit, y=mx+b


gu=length(gh);
%close all; 
 figure
 for i=imshow (sa); % Plots the boundaries of both sides on the imageimshow

  hold on

   plot (col_left,row_left, 'r', 'linewidth', 2);

   plot (col_right,row_right, 'r', 'linewidth', 2);
   
   %plot([qw qt],[qe qr],'r');
   
   plot (gh,bn, 'r', 'linewidth', 2);


 end 

 % Repeat above process with boundaries of breast
w=w_a:w_b;
w=w';

n=numel(w);
c=(y1);
r=repmat(c,n);
r=r(:,2);
 
w_c= min(col_left)
w_d=max(col_left)

w_e= w_c:w_d;
w_e=w_e';

n_i=numel(w_e);
c_i=(y1);
r_i=repmat(c_i,n_i);
r_i=r_i(:,2);

ng=numel(gh);
cv=(y1);
rv=repmat(cv,ng);
rv=rv(:,2);
figure
imshow(sa)

hold on;
for i = 1:length(r)
 plot([w(i),col_right(i)],[r(i),row_right(i)], '-k');
end

for i_i=1:numel(r_i)
  
    plot([w_e(i_i),col_left(i_i)],[r_i(i_i),row_left(i_i)], '-k');
    

end

for i = 1:length(gh)
 plot([gh(i),gh(i)],[rv(i),bn(i)], '-k');
end

hold off;