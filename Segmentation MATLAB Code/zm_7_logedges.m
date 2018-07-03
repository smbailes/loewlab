% Laplacian Using Gaussian Edge Detection
% This code uses the Laplacian of Gaussian method for edge detection. It gives good sides
% and outer lower breast boundaries. It then isolates the lower boundary curves and goes to
% zm_8_curvature.m to find the curvature.

% Last Updated by Zainab Mahmood 8/10/17 at 3:10 PM



%% Part 1: Find 'LoG edges' (laplacian)

% clear
% 
% clc
% close all
% 
% If not using rest of code, uncomment this part so this section can be run
% separately! This part is included in kefset.m and keflarge1
% 
% ptID = input('Enter image name you want to open: ','s'); % Request patient image name
% ptID = strcat(ptID,'.tif'); 
% I = imread(ptID); % open the image, keeping it in 16-bits
% [img_y, img_x] = size(I);
% 
% figure, imshow(I,[]) % to help decide if it should be cropped or not
% title('Image 1: Original Image')

% detect edges using laplacian of gaussian method
% BW = edge(I,'log');
% 
% BW_long = bwareaopen(BW,20);
% BW_long = bwmorph(BW_long,'thicken');
% 
% 
figure;
imshow(I,[]);
title('Log Edges');
magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %magenta has RGB value of 1 0 1
hold on 
h3 = imshow(magenta); 
set(h3, 'AlphaData', BW_long);
hold off


%% Part 2: Overlap final and log edges images
% 
logfin = zeros(img_y,img_x);
for al = 1:img_y
    for bl = 1:img_x
        if newboundaries1(al,bl)==1
            logfin(al,bl)=1;
        end
        if BW_long(al,bl)==1
            logfin(al,bl)=1;
        end
    end
end

figure;
imshow(I,[]);
title('Log and Final');
green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
hold on
h4 = imshow(green);
set(h4, 'AlphaData', logfin);
hold off

% skeleton then connect lines to clean lines up
connected = bwmorph(logfin,'skel',Inf);
connected = bwmorph(connected,'bridge');
connected = bwmorph(connected,'thicken');
connected = bwareaopen(connected,40);
%connected = bwmorph(connected,'clean');

figure;
imshow(I,[]);
title('Log and Final c,s');
green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
hold on
h5 = imshow(green);
set(h5, 'AlphaData', connected);
hold off

% kefpart_7finishconnections

%% Part 3: Graph thicker log and final edges
[dist,dind]=bwdist(BW);
gettem=find(dist<2&dist>0);
if gettem==0
    gettem=1;
end
BW(gettem)=1;

%removes all connected components that have fewer than P=250 pixels 
BW_lo = bwareaopen(BW,250);

figure;
imshow(I,[]);
title('Log and Final Thicker');
cyan = cat(3, zeros(size(I)), ones(size(I)), ones(size(I))); %cyan is 011
hold on
h5 = imshow(cyan);
set(h5, 'AlphaData', BW_lo);


%% Part 4: Split image in half and find connected components

% find connected components
CC_log = bwconncomp(BW_lo);

comp_x = cell(1,CC_log.NumObjects);     % create cell arrays for x and y coordinates for each conn comp
comp_y = cell(1,CC_log.NumObjects);
for n = 1:CC_log.NumObjects
    [comp_y{n}, comp_x{n}] = ind2sub(size(I),CC_log.PixelIdxList{n});               % turn pixel indices into x,y coordinates to graph later
end

% graph components separately
conncompmat = cell(1,CC_log.NumObjects);
for n = 1:CC_log.NumObjects
    conncompmat{n} = zeros(img_y,img_x);
    for a = 1:length(comp_x{n})                                % put each component in a matrix, which is part of a cell array
        conncompmat{n}(comp_y{n}(a),comp_x{n}(a)) = 1;         % each cell is a different component
    end
    
    % graph components
    figure;
    imshow(I,[]);
    mytitle = sprintf('Connected Component # %i',n);
    title(mytitle);
    red = cat(3, ones(size(I)), zeros(size(I)), zeros(size(I))); %red is 100
    hold on
    h7 = imshow(red);
    set(h7, 'AlphaData', conncompmat{n});
    hold off
end

% split image in half
figure;
imshow(I,[]);
hold on;
mid_col = zeros(img_y,img_x);
mid_col(:,img_x/2) = 1;
h6 = imshow(cyan);
set(h6, 'AlphaData', mid_col);
hold off


%% Part 5: Select which component is breast on each side of mid column by first finding closest pixel to midpoint

% BW_lo is the matrix with log edges and center line
mid_ind = find(mid_col==1);
[mid_y, mid_x] = ind2sub(size(I),mid_ind);
mid_coords = horzcat(mid_y,mid_x);
midpt = length(mid_ind)/2;

% make another matrix with log edges and a point in the center of the mid
% line so we can get distances from there and then find the minimum
% BW_mdpt = zeros(img_y, img_x);
% for aa = 1:img_y
%     for bb = 1:img_x
%         if BW_lo(aa,bb)==1
%             BW_mdpt(aa,bb)=1;
%         end
%     end
% end
% BW_mdpt(midpt,mid_coords(1,2))=1;
% 
% figure;
% imshow(I,[]);
% title('Log edges with midpoint');
% red = cat(3, ones(size(I)), zeros(size(I)), zeros(size(I))); %red is 100
% hold on
% h8 = imshow(red);
% set(h8, 'AlphaData', BW_mdpt);
% hold off

% make 2 separate images for left and right breast and 0 out everything
% else

BW_l = zeros(img_y, img_x);     %zero out right breast
for aa = 1:img_y
    for bb = 1:img_x
        if BW_lo(aa,bb)==1
            BW_l(aa,bb)=1;
        end
    end
end
BW_l(:,mid_coords(1,2):end)=0;

BW_r = zeros(img_y, img_x);     %zero out left breast
for aa = 1:img_y
    for bb = 1:img_x
        if BW_lo(aa,bb)==1
            BW_r(aa,bb)=1;
        end
    end
end
BW_r(:,1:mid_coords(1,2))=0;

% find left breast closest pixel
[D_l, IDX_l] = bwdist(BW_l);
[IDX_ly, IDX_lx] = ind2sub(size(I),IDX_l);
IDX_lcoords = zeros(1,2);
IDX_lcoords(1,1) = IDX_ly(midpt,mid_coords(1,2));
IDX_lcoords(1,2) = IDX_lx(midpt,mid_coords(1,2));

% find right breast closest pixel
[D_r, IDX_r] = bwdist(BW_r);
[IDX_ry, IDX_rx] = ind2sub(size(I),IDX_r);
IDX_rcoords = zeros(1,2);
IDX_rcoords(1,1) = IDX_ry(midpt,mid_coords(1,2));
IDX_rcoords(1,2) = IDX_rx(midpt,mid_coords(1,2));

% emphasize pixel found
lpix = zeros(img_y, img_x);
lpix(IDX_lcoords(1),IDX_lcoords(2))=1;
lpix = bwmorph(lpix,'thicken');

rpix = zeros(img_y, img_x);
rpix(IDX_rcoords(1),IDX_rcoords(2))=1;
rpix = bwmorph(rpix,'thicken');

% graph left and right breasts separately and emphasize pixel found
figure;
imshow(I,[]);
title('Log edges left');
red = cat(3, ones(size(I)), zeros(size(I)), zeros(size(I))); %red is 100
hold on
h8 = imshow(red);
set(h8, 'AlphaData', BW_l);
green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
h10 = imshow(green);
set(h10, 'AlphaData', lpix);
hold off

figure;
imshow(I,[]);
title('Log edges right');
red = cat(3, ones(size(I)), zeros(size(I)), zeros(size(I))); %red is 100
hold on
h9 = imshow(red);
set(h9, 'AlphaData', BW_r);
green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
h11 = imshow(green);
set(h11, 'AlphaData', rpix);
hold off

%% Part 6: Find which connected component from conncompmat has that pixel

% find left breast component
for n = 1:CC_log.NumObjects
    if conncompmat{n}(IDX_lcoords(1),IDX_lcoords(2))==1
        foundcomp_l = n;
        fprintf('The left breast component is component # %i \n',foundcomp_l);
        break;
    else
        continue;
    end
end

% find right breast component
for nn = 1:CC_log.NumObjects
    if conncompmat{nn}(IDX_rcoords(1),IDX_rcoords(2))==1
        foundcomp_r = nn;
        fprintf('The right breast component is component # %i \n',foundcomp_r);
        break;
    else
        continue;
    end
end

comp_xnow = comp_x{foundcomp_l};
comp_ynow = comp_y{foundcomp_l};

%% Part 7: Find function for breast component

% fit a fourier curve to breast pixels and exclude outliers
% this code gives the same result as the Curve Fitting GUI

% options = fitoptions('Method','SmoothingSpline','SmoothingParam',0.7);
% f_l = fit(comp_xnow, comp_ynow, 'SmoothingSpline',options)
% figure;
% plot(f_l, comp_xnow, comp_ynow,'.');
% set(gca,'YDir','reverse');
% title('Y Axis Reversed');
% 
% 
% comp_xnowr = comp_x{foundcomp_r};
% comp_ynowr = comp_y{foundcomp_r};
% options = fitoptions('SmoothingSpline');
% f_r = fit(comp_xnowr, comp_ynowr, 'SmoothingSpline',options)
% figure;
% plot(f_r, comp_xnowr, comp_ynowr,'.');
% set(gca,'YDir','reverse');
% title('Y Axis Reversed');

%% Part 8: Calculate curvature

% this gives one curvature value for the whole curve - we want curvature at
% different points so we can track it as a function of position

% x = comp_xnow;
% y = comp_ynow;
% 
% 
% mx = mean(x); my = mean(y);
% X = x - mx; Y = y - my; % Get differences from means
% dx2 = mean(X.^2); dy2 = mean(Y.^2); % Get variances
% t = [X,Y]\(X.^2-dx2+Y.^2-dy2)/2; % Solve least mean squares problem
% a0 = t(1); b0 = t(2); % t is the 2 x 1 solution array [a0;b0]
% r = sqrt(dx2+dy2+a0^2+b0^2); % Calculate the radius
% a = a0 + mx; b = b0 + my; % Locate the circle's center
% curv = 1/r; % Get the curvature

% go on to curvature code
zm_8_curvature;

