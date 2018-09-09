
function [total1 colright colleft] =  connectPixels(total, img_y, img_x, I)

% split image in half
mid_col = zeros(img_y,img_x);
mid_col(:,img_x/2) = 1;

%find coordinates of midline
[midx midy] = find(mid_col==1); 


%find x and y locations of pixels (1s) in logical matrix
[xlocs ylocs] = find(total == 1);
%%[midxx midyy] = find(ylocs == midy(1,1));


[~,c2mid] = min(abs(ylocs-midy(1,1)));
closest = ylocs(c2mid);
midxx = max(find(ylocs==closest));

%% from col 1 to midline col, finds x and y location of lowest pixel

 r = 100;
% 
% %get lowest pixel (row of lower boundary of breast)
% maxx = max(xlocs(1:midxx));
% maxxloc = max(find(xlocs==maxx));
% maxy = ylocs(find(xlocs==maxx)); %finds column of lower boundary of breast (is this necessary??)
% maxy = max(maxy);
% pixx = 0;
% %figure,imshow(total), hold on; plot(maxx,maxy);
% 
% maxx2 = max(xlocs(midxx:end));
% maxxloc2 = max(find(xlocs==maxx2));
% maxy2 = ylocs(find(xlocs==maxx2)); %finds column of lower boundary of breast (is this necessary??)
% maxy2 = max(maxy2);
% pixx = 0;

fprintf('Pick lowest point on the left breast. \n');
figure, imshow(I, [])
[rowleft colleft] = ginput(1);
    
fprintf('Pick lowest point on the right breast. \n');
figure, imshow(I, [])
[rowright colright] = ginput(1);

[~,ml] = min(abs(ylocs-rowleft));
cl = ylocs(ml);
[~,mr] = min(abs(ylocs-rowright));
cr = ylocs(mr);

maxxloc = max(find(ylocs==cr));
maxxloc2 = max(find(ylocs==cl));

%% 
total1 = total; 

%% Section 1

%loop through pixels and find closest pixels 
%while pixx~=maxx %runs until it hits pixel at row of breast lower boundary
%figure, imshow(total), hold on;
%    for xx = 1:maxxloc %runs thru each pixel
     for xx = 1:maxxloc 
         pixx = xlocs(xx,1); %x location of pixel
         pixy = ylocs(xx,1); %y location of pixel
         theta = 0: pi/100 : pi; %angles of lower half circle
         xcirc = r * cos(theta) + pixy; %draws half circle of radius r around pixel (x components)
         ycirc = r * sin(theta) + pixx; %draws half circle of radius r around pixel (y components)
%        plot(xcirc,ycirc), hold on;
        xfound = 0;
        yfound = 0;
        for i = 1:r %loops through first half of points on circle
            bound1x = round(xcirc(i)); %first bound x
                if bound1x<1, bound1x = 1; end
                if bound1x>640, bound1x = 640; end
            bound1y = round(ycirc(i)); %first bound y
                if bound1y<1, bound1y = 1; end
                if bound1y>480, bound1y = 480; end
            bound2x = round(xcirc(length(ycirc)+1-i)); %second bound x
                if bound2x<1, bound2x = 1; end
                if bound2x>640, bound2x = 640; end
            bound2y = round(ycirc(length(ycirc)+1-i)); %second bound y
                if bound2y<1, bound2y = 1; end
                if bound2y>480, bound2y = 480; end
            [yfoundnew] = find(total(bound1y, bound2x:bound1x)==1); %find where there is a pixel on line between first and second bounds 
            [xfoundnew] = ones(1,length(yfoundnew))*(bound1y);
            if yfoundnew 
                xfound = [xfound,xfoundnew];
                yfound = [yfound,yfoundnew];
            end
        end 
        yfound(find(yfound==0)) = [];
        xfound(find(yfound==0)) = [];
        xfound(find(xfound==0)) = [];
        yfound(find(xfound==0)) = [];
        yrem = find(yfound==pixy);
        for j = 1:length(yrem)
            if xfound(yrem)==pixx
                yfound(yrem) = [];
                xfound(yrem) = []; 
            end
        end
        dist = zeros(length(xfound),3);
%       if xfound
            for i = 1:length(xfound) %calculates distance of pixels found in radius from initial pixel
                %if xfound(i)~=0 && yfound(i)~=0 
                    dist(i,1) = sqrt(((yfound(i)-pixy)^2)+((xfound(i)-pixx)^2)); %saves distance to first column
                    dist(i,2) = xfound(i); %saves x location of found pixel to second column
                    dist(i,3) = yfound(i); %saves y location of found pixel to third column
                %end
            end
            [mindistx mindisty] = find(dist(:,1)==min(dist(:,1))); %finds smallest distance and saves location to minddist
            closelocx = max(dist(mindistx,2)); %finds x component of closest found pixel 
            closelocy = max(dist(mindistx,3)); %finds y component of closest found pixel 
            if ~isempty(closelocx) && ~isempty(closelocy)
%                 [xPts,yPts] = bresenham(pixx,pixy,closelocx,closelocy); %finds pixels in between pixel and closest found pixel'
                total1 = linept(total1,pixx,pixy,closelocx,closelocy);
            end
%             yPts(find(yPts==0)) = [];
%             xPts(find(yPts==0)) = [];
%             xPts(find(xPts==0)) = [];
%             yPts(find(xPts==0)) = [];
%             total1(xPts,yPts) = 1; %sets pixels from above to one 
%         end
    end
%    hold off;
%end
%% Section 2
%loop through pixels and find closest pixels 
%while pixx~=maxx %runs until it hits pixel at row of breast lower boundary
%figure, imshow(total), hold on;
    for xx = maxxloc+1:midxx %runs thru each pixel
         pixx = xlocs(xx,1); %x location of pixel
         pixy = ylocs(xx,1); %y location of pixel
         theta = pi: pi/100 : 2*pi; %angles of lower half circle
         xcirc = r * cos(theta) + pixy; %draws half circle of radius r around pixel (x components)
         ycirc = r * sin(theta) + pixx; %draws half circle of radius r around pixel (y components)
%        plot(xcirc,ycirc), hold on;
        xfound = 0;
        yfound = 0;
        for i = 1:r %loops through first half of points on circle
            bound1x = round(xcirc(i)); %first bound x
                if bound1x<1, bound1x = 1; end
                if bound1x>640, bound1x = 640; end
            bound1y = round(ycirc(i)); %first bound y
                if bound1y<1, bound1y = 1; end
                if bound1y>480, bound1y = 480; end
            bound2x = round(xcirc(length(ycirc)+1-i)); %second bound x
                if bound2x<1, bound2x = 1; end
                if bound2x>640, bound2x = 640; end
            bound2y = round(ycirc(length(ycirc)+1-i)); %second bound y
                if bound2y<1, bound2y = 1; end
                if bound2y>480, bound2y = 480; end
 %           [xfoundnew yfoundnew] = find(total(bound1y, bound2x:bound1x)==1); %find where there is a pixel on line between first and second bounds 
            rowsearch = bound1y;
%             [yfoundnew] = find(total(bound1y, bound2x:bound1x)==1); %find where there is a pixel on line between first and second bounds 
%             [xfoundnew] = ones(1,length(yfoundnew))*(bound1y);
%             if yfoundnew 
%                 xfound = [xfound,xfoundnew];
%                 yfound = [yfound,yfoundnew];
%             end
            for colsearch = bound2x:bound1x
                if total(rowsearch, colsearch)==1
                    xfoundnew = rowsearch;
                    yfoundnew = colsearch;
                    xfound = [xfound,xfoundnew];
                    yfound = [yfound,yfoundnew];
                end
            end
%             [xfoundnew] = [];
%             [xfoundnew] = ones(1,length(yfoundnew))*(bound1y);
            %if yfoundnew>=pixy 
                
            %end
        end
        
        yfound(find(yfound==0)) = [];
        xfound(find(yfound==0)) = [];
        xfound(find(xfound==0)) = [];
        yfound(find(xfound==0)) = [];
        
%         xfound(find(yfound<pixy)) = [];
%         yfound(find(yfound<pixy)) = [];

        
%         yrem = find(yfound==pixy);
%         for j = 1:length(yrem)
%             if xfound(yrem)==pixx
%                 yfound(yrem) = [];
%                 xfound(yrem) = []; 
%             end
%         end
        dist = zeros(length(xfound),3);
%       if xfound
            for i = 1:length(xfound) %calculates distance of pixels found in radius from initial pixel
                %if xfound(i)~=0 && yfound(i)~=0 
                    dist(i,1) = sqrt(((yfound(i)-pixy)^2)+((xfound(i)-pixx)^2)); %saves distance to first column
                    dist(i,2) = xfound(i); %saves x location of found pixel to second column
                    dist(i,3) = yfound(i); %saves y location of found pixel to third column
                %end
            end
            [mindistx mindisty] = find(dist(:,1)==min(dist(:,1))); %finds smallest distance and saves location to minddist
            closelocx = max(dist(mindistx,2)); %finds x component of closest found pixel 
            closelocy = max(dist(mindistx,3)); %finds y component of closest found pixel 
            if ~isempty(closelocx) && ~isempty(closelocy)
%                 [xPts,yPts] = bresenham(pixx,pixy,closelocx,closelocy); %finds pixels in between pixel and closest found pixel'
                  total1 = linept(total1,pixx,pixy,closelocx,closelocy);
            end
%             yPts(find(yPts==0)) = [];
%             xPts(find(yPts==0)) = [];
%             xPts(find(xPts==0)) = [];
%             yPts(find(xPts==0)) = [];
%             total1(xPts,yPts) = 1; %sets pixels from above to one 
%         end
%         closelocx
%         closelocy
    end
%    hold off;
%end

%% Section 3
%loop through pixels and find closest pixels 
%while pixx~=maxx %runs until it hits pixel at row of breast lower boundary
%figure, imshow(total), hold on;
    for xx = midxx:maxxloc2 %runs thru each pixel
         pixx = xlocs(xx,1); %x location of pixel
         pixy = ylocs(xx,1); %y location of pixel
         theta = 0: pi/100 : pi; %angles of lower half circle
         xcirc = r * cos(theta) + pixy; %draws half circle of radius r around pixel (x components)
         ycirc = r * sin(theta) + pixx; %draws half circle of radius r around pixel (y components)
%        plot(xcirc,ycirc), hold on;
        xfound = 0;
        yfound = 0;
        for i = 1:r %loops through first half of points on circle
            bound1x = round(xcirc(i)); %first bound x
                if bound1x<1, bound1x = 1; end
                if bound1x>640, bound1x = 640; end
            bound1y = round(ycirc(i)); %first bound y
                if bound1y<1, bound1y = 1; end
                if bound1y>480, bound1y = 480; end
            bound2x = round(xcirc(length(ycirc)+1-i)); %second bound x
                if bound2x<1, bound2x = 1; end
                if bound2x>640, bound2x = 640; end
            bound2y = round(ycirc(length(ycirc)+1-i)); %second bound y
                if bound2y<1, bound2y = 1; end
                if bound2y>480, bound2y = 480; end
 %           [xfoundnew yfoundnew] = find(total(bound1y, bound2x:bound1x)==1); %find where there is a pixel on line between first and second bounds 
            rowsearch = bound1y;
            for colsearch = bound2x:bound1x
                if total(rowsearch, colsearch)==1
                    xfoundnew = rowsearch;
                    yfoundnew = colsearch;
                    xfound = [xfound,xfoundnew];
                    yfound = [yfound,yfoundnew];
                end
            end
%             [xfoundnew] = [];
%             [xfoundnew] = ones(1,length(yfoundnew))*(bound1y);
            %if yfoundnew>=pixy 
                
            %end
        end
        
        yfound(find(yfound==0)) = [];
        xfound(find(yfound==0)) = [];
        xfound(find(xfound==0)) = [];
        yfound(find(xfound==0)) = [];
        
%         xfound
%         yfound
        
%         xfound(find(yfound<pixy)) = [];
%         yfound(find(yfound<pixy)) = [];
% 
%         
%         yrem = find(yfound==pixy);
%         for j = 1:length(yrem)
%             if xfound(yrem)==pixx
%                 yfound(yrem) = [];
%                 xfound(yrem) = []; 
%             end
%         end
        dist = zeros(length(xfound),3);
%       if xfound
            b = 1;
            for i = 1:length(xfound) %calculates distance of pixels found in radius from initial pixel
                %if xfound(i)~=0 && yfound(i)~=0 
                        dist(i,1) = sqrt(((yfound(i)-pixy)^2)+((xfound(i)-pixx)^2)); 
                        dist(i,2) = xfound(i); %saves x location of found pixel to second column
                        dist(i,3) = yfound(i); %saves y location of found pixel to third column
                %end
            end
%              dist(find(dist(:,1)==0)) = [];
%              dist(find(dist(:,2)==0)) = [];
%              dist(find(dist(:,2)==0)) = [];
            [rowrem colrem] = find((dist(:,1)==0));
            dist(rowrem,:)=[];
            [mindistx mindisty] = find(dist(:,1)==min(dist(:,1))); %finds smallest distance and saves location to minddist
            closelocx = max(dist(mindistx,2)); %finds x component of closest found pixel 
            closelocy = max(dist(mindistx,3)); %finds y component of closest found pixel 
            if ~isempty(closelocx) && ~isempty(closelocy)
%                 [xPts,yPts] = bresenham(pixx,pixy,closelocx,closelocy); %finds pixels in between pixel and closest found pixel'
                  total1 = linept(total1,pixx,pixy,closelocx,closelocy);
            end
%             yPts(find(yPts==0)) = [];
%             xPts(find(yPts==0)) = [];
%             xPts(find(xPts==0)) = [];
%             yPts(find(xPts==0)) = [];
%             total1(xPts,yPts) = 1; %sets pixels from above to one 
%         end
%     figure, imshow(total1)
%     pixx
%     pixy
%     closelocx
%     closelocy
%     dist
%     break;
    end
%    hold off;
%end

%% Section 4

%loop through pixels and find closest pixels 
%while pixx~=maxx %runs until it hits pixel at row of breast lower boundary
%figure, imshow(total), hold on;
len = length(xlocs);
    for xx = maxxloc2+1:len %runs thru each pixel
         pixx = xlocs(xx,1); %x location of pixel
         pixy = ylocs(xx,1); %y location of pixel
         theta = pi: pi/100 : 2*pi; %angles of lower half circle
         xcirc = r * cos(theta) + pixy; %draws half circle of radius r around pixel (x components)
         ycirc = r * sin(theta) + pixx; %draws half circle of radius r around pixel (y components)
%        plot(xcirc,ycirc), hold on;
        xfound = 0;
        yfound = 0;
        for i = 1:r %loops through first half of points on circle
            bound1x = round(xcirc(i)); %first bound x
                if bound1x<1, bound1x = 1; end
                if bound1x>640, bound1x = 640; end
            bound1y = round(ycirc(i)); %first bound y
                if bound1y<1, bound1y = 1; end
                if bound1y>480, bound1y = 480; end
            bound2x = round(xcirc(length(ycirc)+1-i)); %second bound x
                if bound2x<1, bound2x = 1; end
                if bound2x>640, bound2x = 640; end
            bound2y = round(ycirc(length(ycirc)+1-i)); %second bound y
                if bound2y<1, bound2y = 1; end
                if bound2y>480, bound2y = 480; end
 %           [xfoundnew yfoundnew] = find(total(bound1y, bound2x:bound1x)==1); %find where there is a pixel on line between first and second bounds 
            rowsearch = bound1y;
            for colsearch = bound2x:bound1x
                if total(rowsearch, colsearch)==1
                    xfoundnew = rowsearch;
                    yfoundnew = colsearch;
                    xfound = [xfound,xfoundnew];
                    yfound = [yfound,yfoundnew];
                end
            end
%             [xfoundnew] = [];
%             [xfoundnew] = ones(1,length(yfoundnew))*(bound1y);
            %if yfoundnew>=pixy 
                
            %end
        end
        
        yfound(find(yfound==0)) = [];
        xfound(find(yfound==0)) = [];
        xfound(find(xfound==0)) = [];
        yfound(find(xfound==0)) = [];
        
%         xfound(find(yfound<pixy)) = [];
%         yfound(find(yfound<pixy)) = [];
% 
%         
%         yrem = find(yfound==pixy);
%         for j = 1:length(yrem)
%             if xfound(yrem)==pixx
%                 yfound(yrem) = [];
%                 xfound(yrem) = []; 
%             end
%         end
        dist = zeros(length(xfound),3);
%       if xfound
            for i = 1:length(xfound) %calculates distance of pixels found in radius from initial pixel
                %if xfound(i)~=0 && yfound(i)~=0 
                    dist(i,1) = sqrt(((yfound(i)-pixy)^2)+((xfound(i)-pixx)^2)); %saves distance to first column
                    dist(i,2) = xfound(i); %saves x location of found pixel to second column
                    dist(i,3) = yfound(i); %saves y location of found pixel to third column
                %end
            end
            [mindistx mindisty] = find(dist(:,1)==min(dist(:,1))); %finds smallest distance and saves location to minddist
            closelocx = max(dist(mindistx,2)); %finds x component of closest found pixel 
            closelocy = max(dist(mindistx,3)); %finds y component of closest found pixel 
            if ~isempty(closelocx) && ~isempty(closelocy)
%                 [xPts,yPts] = bresenham(pixx,pixy,closelocx,closelocy); %finds pixels in between pixel and closest found pixel'
                total1 = linept(total1,pixx,pixy,closelocx,closelocy);
            end
%    hold off;
    end
%% Images

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
end
%% 
% function [x,y]=bresenham(x1,y1,x2,y2)
% % Line to pixel approximation algorithm (Bresenham's line algorithm)
% % Credit: Aaron WetzlerAll (2010)
% x1=round(x1); x2=round(x2);
% y1=round(y1); y2=round(y2);
% dx=abs(x2-x1);dy=abs(y2-y1);
% steep=abs(dy)>abs(dx);
% if steep t=dx;dx=dy;dy=t; end
% 
% if dy==0 
%     q=zeros([dx+1,1]);
% else
%     q=[0;diff(mod((floor(dx/2):-dy:-dy*dx+floor(dx/2))',dx))>=0];
% end
% 
% if steep
%     if y1<=y2 y=(y1:y2)'; else y=(y1:-1:y2)'; end
%     if x1<=x2 x=x1+cumsum(q);else x=x1-cumsum(q); end
% else
%     if x1<=x2 x=(x1:x2)'; else x=(x1:-1:x2)'; end
%     if y1<=y2 y=y1+cumsum(q);else y=y1-cumsum(q); end
% end
% end
%% 
function result=linept(matrice, X1, Y1, X2, Y2)
% Connect two pixels in a matrice with 1
%
% Command line
% ------------
% result=linept(matrice, X1, Y1, X2, Y2)
%   matrice : matrice where I'll write
%   (X1, Y1), (X2, Y2) : points to connect
%   result : matrix + the line
%
% Note
% ----
%   matrice can contents anything
%   (X1, Y1), (X2, Y2) can be out of the matrice
%
% Example
% -------
% a = linept(zeros(5, 10), 2, 2, 3, 9)
% a =
% 
%      0     0     0     0     0     0     0     0     0     0
%      0     1     1     1     1     0     0     0     0     0
%      0     0     0     0     0     1     1     1     1     0
%      0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0
%
% Georges Cubas 20/11/03
% georges.c@netcourrier.com
% Version 1.0

result = matrice;
for x=max(1, X1):sign(X2 - X1):max(1, X2)
    y = round(f(x, X1, Y1, X2, Y2));
    if y > 0
        result(x, y) = 1;
    end
end
for y=max(1, Y1):sign(Y2 - Y1):max(1, Y2)
    x = round(f2(y, X1, Y1, X2, Y2));
    if x > 0
        result(x, y) = 1;
    end
end

function y=f(x, X1, Y1, X2, Y2)
a = (Y2 - Y1)/(X2 - X1);
b = Y1 - X1 * a;
y = a * x + b;
end

function x=f2(y, X1, Y1, X2, Y2)
if X1==X2
    x = X1;
else
	a = (Y2 - Y1)/(X2 - X1);
	b = Y1 - X1 * a;
	x = (y - b)/a;
end
end
end
