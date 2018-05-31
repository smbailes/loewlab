% Curvature calculation and continuation

% This code finds the curvature using rate of turn and replicates it to continue the
% curve on the inner sides of the breasts.

% Last Updated by Zainab Mahmood 8/10/17 at 3:10PM

%% Part 1: Make boundary one pixel thick to find curvature

% remove points that decrease
for r = 2:length(comp_ynow)-1
    if comp_ynow(r-1)<comp_ynow(r) && comp_ynow(r+1)>comp_ynow(r)
        continue;
    else
        comp_ynow(r)=0;
    end
end

findzero = find(comp_ynow==0);
comp_ynow(findzero)=[];
comp_xnow(findzero)=[];

% need to get rid of multiple pixels with the same y
% this just puts pixels of the component into a matrix to use later
tempmat = zeros(img_y, img_x);
for a = 1:length(comp_xnow)
    tempmat(comp_ynow(a),comp_xnow(a))=2;
end


comp_ynow_thin = bwmorph(tempmat,'skel');           % run skel twice to thin boundary even more
comp_ynow_thinn = bwmorph(comp_ynow_thin,'skel',Inf);
comp_ynow_thin = double(comp_ynow_thin);        % cast logical to double
comp_ynow_thinn = double(comp_ynow_thinn);
for cnt = 1:img_y                     % change 1's to a large number that's easier to see in the matrix
    for cnt2 = 1:img_x
        if comp_ynow_thin(cnt,cnt2)==1
            comp_ynow_thin(cnt,cnt2)=2^16;
        end
        if comp_ynow_thinn(cnt,cnt2)==1
            comp_ynow_thinn(cnt,cnt2)=2^16;
        end
    end
end

% make boundary one pixel thick
thinmat = zeros(img_y, img_x);
temp = zeros(1, img_x);
for count = 1:img_x                     % for all columns
    temp(count) = sum(comp_ynow_thinn(:,count));    % get sum of each column
    if temp(count)~=0                               % if sum isn't zero
        loc = find(comp_ynow_thinn(:,count));       % find indices of which rows aren't zero
        if loc(end)-loc(1)>1                        % if pixels aren't next to each other, keep first and last
            thinmat(loc(1),count) = 2^16;
            thinmat(loc(end),count) = 2^16;
        else                                        % otherwise only keep bottom pixel
            thinmat(loc(end),count) = 2^16;
        end
    else
        continue;
    end
end

% get x,y coordinates
findthin = find(thinmat>0);
[thiny, thinx] = ind2sub(size(I),findthin);
startpt = input('Enter starting point for contour: ');
direction = input('Enter cw or ccw: ','s');
[xout, yout] = points2contour(thinx,thiny,startpt,direction);
figure;
comet(xout,yout);
title('Closed contour');
set(gca,'YDir','reverse');

% put x,y coordinates into matrix
thin_bound = zeros(img_y, img_x);
for cnt = 1:length(xout)
    thin_bound(yout(cnt),xout(cnt))=2^16;
end

figure, imshow(I,[]), title('Thin Boundary')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', thin_bound);

% try to make it horizontally one pixel thin
xout2 = xout;
yout2 = yout;
for out = 2:length(xout)
    if yout(out)==yout(out-1)
        yout2(out)=0;
        xout2(out)=0;
    end
end

% take out 0'ed pixels
oo = find(yout2~=0);
yout2 = yout2(oo);
oo2 = find(xout2~=0);
xout2 = xout2(oo2);

% put x,y coordinates into matrix
thin_bound2 = zeros(img_y, img_x);
for cnt = 1:length(xout2)
    thin_bound2(yout2(cnt),xout2(cnt))=2^16;
end

figure, imshow(I,[]), title('Thin Boundary2')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', thin_bound2);

% % smoothing filter
% fs = 640;
% Nyf = fs/2;
% 
% W1 = 40;
% W2 = 50;
% [ind1,W] = buttord(W1/Nyf,W2/Nyf,3,10);
% [indb, inda] = butter(ind1,W);
% filtered = filtfilt(indb,inda,comp_ynow);
% figure;
% subplot(1,2,1);
% plot(comp_xnow,comp_ynow);
% hold on;
% plot(comp_xnow,comp_ynow,'o');
% subplot(1,2,2);
% plot(comp_xnow,filtered);
% 
% % gaussian filter - didn't really do anything
% gausfilt = imgaussfilt(comp_ynow_thinn);
% gausfilt = bwmorph(gausfilt,'skel',Inf);
% gausfilt = double(gausfilt);
% for cnt = 1:img_y                     % change 1's to a large number that's easier to see in the matrix
%     for cnt2 = 1:img_x
%         if gausfilt(cnt,cnt2)==1
%             gausfilt(cnt,cnt2)=2^16;
%         end
%     end
% end

% find curvature using curvature function from FileExchange
% vertices = horzcat(xout',yout');
% curvature = LineCurvature2D(vertices);
% figure;
% % subplot(1,2,1);
% plot(curvature);

%% Part 2: Calculate curvature using rate of turn

yout2 = 480-yout2;
thinpts = horzcat(xout2',yout2');
% find slope btwn next and previous pt and its magnitude
m = zeros(1,length(thinpts));
m_mag = zeros(1,length(thinpts));
for cntr = 2:length(thinpts)-1          % use on all points except for first and last
    y2 = thinpts(cntr-1,2);             % (x2,y2) is next pt, (x1,y1) is previous pt
    y1 = thinpts(cntr+1,2);
    x2 = thinpts(cntr-1,1);
    x1 = thinpts(cntr+1,1);
    m(cntr) = (y2-y1)/(x2-x1);          % calculate slope
    m_mag(cntr) = sqrt((x2-x1)^2+(y2-y1)^2);    % calculate magnitude of slope
end

figure;
plot(m);
title('Slope m');

% find unit tangent vector
t = zeros(1,length(thinpts));
for val = 2:length(thinpts)-1
    t(val) = m(val)./m_mag(val);
end

figure;
plot(t);
title('Unit Tangent Vector t(t)');

% find psi
psi = zeros(1,length(thinpts));
for value = 2:length(thinpts)-1
    y2 = thinpts(value-1,2);             % (x2,y2) is next pt, (x1,y1) is previous pt
    y1 = thinpts(value+1,2);
    x2 = thinpts(value-1,1);
    x1 = thinpts(value+1,1);
    psi(value) = atan((abs(y2-y1))/(abs(x2-x1)));       % calculate slope angle
end

% find change in psi
dpsi = diff(psi);

% find distance (length) btwn next and previous pt
len = zeros(1,length(thinpts));
for cntr2 = 2:length(thinpts)-1
    y2 = thinpts(cntr2-1,2);             % (x2,y2) is next pt, (x1,y1) is previous pt
    y1 = thinpts(cntr2+1,2);
    x2 = thinpts(cntr2-1,1);
    x1 = thinpts(cntr2+1,1);
    len(cntr2) = sqrt((x2-x1)^2 + (y2-y1)^2);       % calculate length
end

figure;
plot(len);
title('Length');

% find change in lengths
dlen = diff(len);

% find kappa (curvature)
kappa = zeros(1,length(thinpts));
for kcount = 2:length(thinpts)-1
    if dlen(kcount)~=0
        kappa(kcount) = abs(dpsi(kcount))/abs(dlen(kcount));
    else
        kappa(kcount) = 0;
    end
end
kk = find(kappa~=0);
kappa2 = kappa(kk);
figure;
plot(kappa2);
title('Curvature');
xlabel('Pixel Position');
ylabel('Curvature');

%% Part 3: Continue curvature

[min_k,min_kind] = min(kappa2);
kappa_new = kappa2(1:min_kind);
kappa_new = [kappa_new fliplr(kappa_new)];

figure;
plot(kappa_new);
title('Curvature continued');
xlabel('Pixel Position');
ylabel('Curvature');

zm_85_topline;

