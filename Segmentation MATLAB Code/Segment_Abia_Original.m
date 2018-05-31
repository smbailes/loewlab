% Segmentation Algorithm
% Author : Abia Khan
% Last revised: 12/13/2016

% The segmentation method is below. 

%% Step 1: Preprocessing: Input Image, Clear background, Revised egde detection algorithm
clear;
close all;
clc;
% Read Image
% sa=imread('/Users/abiakhan/Downloads/Images Loew/0000 (1).tif');
% sb=imread('/Users/abiakhan/Downloads/Images Loew/0001 (1).tif');
% sc=imread('/Users/abiakhan/Downloads/Images Loew/0002 (1).tif');
% sd=imread('/Users/abiakhan/Downloads/Images Loew/0008 (1).tif');
sa=imread_irfanview('1799 - P5.tif');
% Crop Image
sa=imcrop(sa);

%Clear Background
% thresholdvalue = graythresh(sa)*2^16; %Clearing the background : Otsu's method
% 
% [a,b]=size(sa);
% 
% for k = 1:a
% 
%     for u = 1:b
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

imshow(sa)


edge=edge(sa,'canny','vertical'); % Edge Detection, just using the vertical derivative
figure,imshow(edge);
bwImgs = bwareaopen(edge,150); % removes 200 pixes with a connectivity of 8 in the image. This gets rid of almost all of the vessels/spurious lines in the breast.
figure,imshow(bwImgs);

bounds = bwboundaries(bwImgs); % Computes the boundaries in the image


%% Step 2: Display right and left breast boundaries

boundary = bounds{2}; %gives the properties of the second boundary, which is the right breast

plot(boundary(:,2), boundary(:,1),'r'); 

x=boundary(:,1); % rows of right breast boundary

y=boundary(:,2); % columns of right breast boundary

boundaryl = bounds{3};

a=boundaryl(:,1); % rows of left breast boundary

b=boundaryl(:,2); % columns of left breast boundary


 %% Perform Segmentation
 
close all; imshow(sa)
y1=get(gca,'ylim'); % the y limit of the image
close all;
w_a= min(y); %Lowest y coordinate for right boundary
w_b=max(y); %Highest y coordinate for right boundary

qw = max(y);
qe = min(x);

qr = min(a);
qt = min(b);


y_i= [qe qr]; % Create coordinates
x_i= [qw qt];  % Create coordinates


gh = [qw:qt]; %Create columns
gh=gh';
fg=[qr:qe];
fg=fg';
 

f = polyfit(x_i,y_i,1);%Draw best fit line between the region between the boundaries
f=(f(1));

bv = (qr-(f*qt)); % Perform the slope
vb= (f*qw)+bv;

bn = (f*gh)+bv; %Line of best fit, y=mx+b


gu=length(gh);
close all;
 for i=imshow (sa); % Plots the boundaries of both sides on the imageimshow

  hold on
  
   plot (b,a, 'r', 'linewidth', 2);

   plot (y,x, 'r', 'linewidth', 2);
   
   %plot([qw qt],[qe qr],'r');
   
   plot (gh,bn, 'r', 'linewidth', 2);


 end 

 % Repeat above process with boundaries of breast
w=w_a:w_b;
w=w'

n=numel(w);
c=(y1)
r=repmat(c,n);
r=r(:,2)
 
w_c= min(b)
w_d=max(b)

w_e= w_c:w_d;
w_e=w_e'

n_i=numel(w_e);
c_i=(y1)
r_i=repmat(c_i,n_i);
r_i=r_i(:,2)

ng=numel(gh);
cv=(y1)
rv=repmat(cv,ng);
rv=rv(:,2)
imshow(sa)

hold on;
for i = 1:length(r)
 plot([w(i),y(i)],[r(i),x(i)], '-k');
end

for i_i=1:numel(r_i)
  
    plot([w_e(i_i),b(i_i)],[r_i(i_i),a(i_i)], '-k');
    

end

for i = 1:length(gh)
 plot([gh(i),gh(i)],[rv(i),bn(i)], '-k');
end

hold off;






