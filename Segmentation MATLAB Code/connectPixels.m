%% find midline

% split image in half
mid_col = zeros(img_y,img_x);
mid_col(:,img_x/2) = 1;

%find coordinates of midline
[midx midy] = find(mid_col==1); 
midlinez(:,1) = midx;
midlinez(:,2) = midy;

%find x and y locations of pixels (1s) in logical matrix
[xlocs ylocs] = find(total == 1);
[midxx midyy] = find(ylocs == midlinez(1,2));


%% from col 1 to midline col, finds x and y location of lowest pixel

r = 50;

%get lowest pixel (row of lower boundary of breast)
maxx = max(xlocs);
maxxloc = max(find(xlocs==maxx));
maxy = ylocs(find(xlocs==maxx)); %finds column of lower boundary of breast (is this necessary??)
maxy = max(maxy);
pixx = 0;
%% 

total1 = total; 
%% 

%loop through pixels and find closest pixels 
%while pixx~=maxx %runs until it hits pixel at row of breast lower boundary
figure, imshow(total), hold on;
    for xx = 1:maxxloc %runs thru each pixel
         pixx = xlocs(xx,1); %x location of pixel
         pixy = ylocs(xx,1); %y location of pixel
         theta = 0: 0.01 : pi; %angles of lower half circle
         xcirc = r * cos(theta) + pixy; %draws half circle of radius r around pixel (x components)
         ycirc = r * sin(theta) + pixx; %draws half circle of radius r around pixel (y components)
         plot(xcirc,ycirc), hold on;
         for i = 1:r %loops through first half of points on circle
            bound1x = round(xcirc(i)); %first bound x
                if bound1x<1, bound1x = 1; end
                if bound1x>640, bound1x = 640; end
            bound1y = round(ycirc(i)); %first bound y
                if bound1y<1, bound1y = 1; end
                if bound1y>480, bound1y = 480; end
            bound2x = round(xcirc(length(ycirc)-i)); %second bound x
                if bound2x<1, bound2x = 1; end
                if bound2x>640, bound2x = 640; end
            bound2y = round(ycirc(length(ycirc)-i)); %second bound y
                if bound2y<1, bound2y = 1; end
                if bound2y>480, bound2y = 480; end
            [xfound yfound] = find(total(bound1y, bound1x:bound2x)==1); %find where there is a pixel on line between first and second bounds 
         end 
        dist = zeros(length(xfound),3);
        if xfound~=0
            for i = 1:length(xfound) %calculates distance of pixels found in radius from initial pixel
                dist(i,1) = sqrt(((yfound(i)-pixy)^2)+((xfound(i)-pixx)^2)); %saves distance to first column
                dist(i,2) = xfound(i); %saves x location of found pixel to second column
                dist(i,3) = yfound(i); %saves y location of found pixel to third column
            end
            [mindistx mindisty] = find(dist(:,1)==min(dist(:,1))); %finds smallest distance and saves location to minddist
            closelocx = dist(mindistx,2); %finds x component of closest found pixel 
            closelocy = dist(mindistx,3); %finds y component of closest found pixel 
            [xPts,yPts] = bresenham(pixx,pixy,closelocx,closelocy); %finds pixels in between pixel and closest found pixel
            total1(xPts,yPts) = 1; %sets pixels from above to one 
        end
    end
    hold off;
%end
%% 
figure, imshow(I,[]), title('Total Original')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', total)

figure, imshow(I,[]), title('Connect Pixels')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', total1)
    %from first pixel to pixel at x and y location (from above), find closest pixel withing certian radius AND below y comp of that pixel
    %from x and y loc (from above) to midline, finds closest pixel within certain radius AND above y comp of that pixel
    
%% from midline col to last col, finds x and y location of lowest pixel
    %from midline to pixel at x and y location (from above), find closest pixel withing certian radius AND below y comp of that pixel
    %from x and y loc (from above) to last pixel, finds closest pixel within certain radius AND above y comp of that pixel
%% 
function [x,y]=bresenham(x1,y1,x2,y2)
% Line to pixel approximation algorithm (Bresenham's line algorithm)
% Credit: Aaron WetzlerAll (2010)
x1=round(x1); x2=round(x2);
y1=round(y1); y2=round(y2);
dx=abs(x2-x1);dy=abs(y2-y1);
steep=abs(dy)>abs(dx);
if steep t=dx;dx=dy;dy=t; end

if dy==0 
    q=zeros([dx+1,1]);
else
    q=[0;diff(mod((floor(dx/2):-dy:-dy*dx+floor(dx/2))',dx))>=0];
end

if steep
    if y1<=y2 y=(y1:y2)'; else y=(y1:-1:y2)'; end
    if x1<=x2 x=x1+cumsum(q);else x=x1-cumsum(q); end
else
    if x1<=x2 x=(x1:x2)'; else x=(x1:-1:x2)'; end
    if y1<=y2 y=y1+cumsum(q);else y=y1-cumsum(q); end
end
end
    
    