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

r = 25;

%get lowest pixel
maxx = max(xlocs); 
maxy = ylocs(find(xlocs==maxx)); %???

pixx = 0;
%% 

%loop through pixels and find closest pixels 
while pixx<maxx
    for xx = 1:length(xlocs)
         pixx = xlocs(xx,1);
         pixy = ylocs(xx,1); 
         theta = pi: 0.01 : 2*pi;
         xcirc = r * cos(theta) + pixx;
         ycirc = r * sin(theta) + pixy;
         for i = 1:r
            bound1x = round(xcirc(i));
            bound1y = round(ycirc(i));
            bound2x = round(xcirc(length(ycirc)-i)); 
            bound2y = round(ycirc(length(ycirc)-i));
            [xfound yfound] = find(I(bound1x:bound2x, bound1y:bound2y)==1);
         end 
        if xfound~=0
            for i = 1:length(xfound)
                dist(i) = sqrt(((yfound-pixy)^2)+((xfound-pixx)^2));
            end
            [closelocx closelocy] = find(dist==min(dist));
            [xPts,yPts] = bresenham(pixx,pixy,closelocx,closelocy);
            total(xPts,yPts) = 1; 
        end
    end
end
%% 

figure, imshow(I,[]), title('Connect Pixels')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', total)
    %from first pixel to pixel at x and y location (from above), find closest pixel withing certian radius AND below y comp of that pixel
    %from x and y loc (from above) to midline, finds closest pixel within certain radius AND above y comp of that pixel
    
%% from midline col to last col, finds x and y location of lowest pixel
    %from midline to pixel at x and y location (from above), find closest pixel withing certian radius AND below y comp of that pixel
    %from x and y loc (from above) to last pixel, finds closest pixel within certain radius AND above y comp of that pixel
%% 
function [x,y]=bresenham(x1,y1,x2,y2)
% Line to pixel approximation algorithm (Bresenham's line algorithm)
% Credit: Aaron WetzlerAll (2010)
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
    
    