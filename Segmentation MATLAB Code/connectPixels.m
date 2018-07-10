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

%get lowest pixel
minx = min(xlocs); 
miny = ylocs(minx); %???

%loop through pixels and find closest pixels 
for yy = 1:midlinez(1,2)
     for xx = 1:xlocs(midyy)
         pixx = xlocs(xx,yy);
         pixy = ylocs(xx,yy); 
         if (pix
             
         
     end
end
    %from first pixel to pixel at x and y location (from above), find closest pixel withing certian radius AND below y comp of that pixel
    %from x and y loc (from above) to midline, finds closest pixel within certain radius AND above y comp of that pixel
    
%% from midline col to last col, finds x and y location of lowest pixel
    %from midline to pixel at x and y location (from above), find closest pixel withing certian radius AND below y comp of that pixel
    %from x and y loc (from above) to last pixel, finds closest pixel within certain radius AND above y comp of that pixel
    
    