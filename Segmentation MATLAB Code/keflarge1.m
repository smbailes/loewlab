% keflarge1
%%  Part 1, Large Points
largepoints=zeros(img_y,img_x);


BW = edge(I,'log');

%removes all connected components that have fewer than P=30 pixels 
BW_long = bwareaopen(BW,30); 
%thickens objects by adding pixels to the exterior of objects
BW_long = bwmorph(BW_long,'thicken');

%implementing combined point systems
for aa = 1:img_y
    for bb = 1:img_x
        if edgecanny2(aa,bb)==1%add point if canny edge
            largepoints(aa,bb)=largepoints(aa,bb)+1; 
        end
        if ellipses(aa,bb)==1 %add points if included in ellipse
            largepoints(aa,bb)=largepoints(aa,bb)+2;
        end
        if newI(aa,bb)~=0 %add point for hot pixel
            largepoints(aa,bb)=largepoints(aa,bb)+1;
        end
        if BW_long(aa,bb)~=0 %add points for LoG edge detection  
            largepoints(aa,bb)=largepoints(aa,bb)+2;
        end
    end
end

figure,imshow(largepoints,[])
title('Combined Point Systems')

overlayedpoints = zeros(img_y,img_x); 
for cc=1:img_y 
    for dd=1:img_x
        if largepoints(cc,dd)>2 
            overlayedpoints(cc,dd)=1;
        end
    end
end

[e,f]=size(ellipses);
bottoms=zeros(e,f);
for fir=1:f
    y=find(overlayedpoints(:,fir)==1);
    if ~isempty(y)
    bottoms(y(end),fir)=1;
    end
end

% zainab edit
bottoms2 = double(bottoms);        % cast logical to double
for cnt = 1:img_y                     % change 1's to a large number that's easier to see in the matrix
    for cnt2 = 1:img_x
        if bottoms(cnt,cnt2)==1
            bottoms2(cnt,cnt2)=2^16;
        end
    end
end

%make boundary one pixel thick
bottoms3 = zeros(img_y, img_x);
temp = zeros(1, img_x);
for count = 1:img_x                     % for all columns
    temp(count) = sum(bottoms2(:,count));    % get sum of each column
    if temp(count)~=0                               % if sum isn't zero
        loc = find(bottoms2(:,count));       % find indices of which rows aren't zero
        if loc(end)-loc(1)>1                        % if pixels aren't next to each other, keep first and last
            bottoms3(loc(1),count) = 2^16;
            bottoms3(loc(end),count) = 2^16;
        else                                        % otherwise only keep bottom pixel
            bottoms3(loc(end),count) = 2^16;
        end
    else
        continue;
    end
end

findthin = find(bottoms2>0);
[thiny, thinx] = ind2sub(size(I),findthin);
[xout, yout] = points2contour(thinx,thiny,1,'ccw');

bottoms3 = zeros(img_y, img_x);
for cnt = 1:length(xout)
    bottoms3(yout(cnt),xout(cnt))=2^16;
end

% % 
uppers=zeros(e,f);
for fir=1:f
    y=find(overlayedpoints(:,fir)==1);
    if ~isempty(y)
    uppers(y(1),fir)=1;
    end
end

findthin2 = find(uppers>0);
[thiny2, thinx2] = ind2sub(size(I),findthin2);
[xout2, yout2] = points2contour(thinx2,thiny2,1,'ccw');

uppers2 = zeros(img_y, img_x);
for cnt = 1:length(xout2)
    uppers2(yout2(cnt),xout2(cnt))=2^16;
end

figure, imshow(I,[]), title('Overlayed')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', overlayedpoints)

figure, imshow(I,[]), title('Lowers')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', bottoms3)

figure, imshow(I,[]), title('Uppers')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off
set(displ, 'AlphaData', uppers2)

%total=bwmorph(total,'clean'); %removes individual 1's surrounded by 0's

side=input('Are the breasts lower or higher? [h/l]: ','s');
  if side=='l'
    total = bottoms3;
  else
    total = uppers2;
  end
% 
imwrite(total, 'Total_0000_P8.tif');
% 
figure, imshow(I,[]), title('Best Fit')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', total)

   figure, imshow(I,[]), title('Best Fit')
   % blue on top on figure
   blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
   hold on 
   displ = imshow(blue); 
   hold off 
   set(displ, 'AlphaData', total)
  
%% Part 2, Clean

gettit=bwmorph(total,'close');

gett=bwmorph(gettit,'bridge');

sg = strel('disk',4); %Create a Morphological structuring element, you change the shape used and diameter
gett= imclose(gett,sf);

gett=bwmorph(gett,'clean');

gett=bwmorph(gett,'bridge');


figure, imshow(I,[]), title('Cleaned')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on
displ = imshow(blue);
hold off
set(displ, 'AlphaData', gett)

%% Part 3, Connect

newboundaries = connectDots(gett,50);

% CC = bwconncomp(newboundaries);
% newboundaries = newboundaries;
%  
% for n = 1:CC.NumObjects - 1 %the number of lines in the middle region of the patient
% %     Store all row and col values of component n and the component after in
% %     x1,y1, x2, y2
%     [x1, y1] = ind2sub(size(newboundaries),CC.PixelIdxList{n});
%     [x2, y2] = ind2sub(size(newboundaries),CC.PixelIdxList{n+1});
%     
%     [yy1, ind] = max(y1); %find the max col in component n
%     xx1 = x1(ind); % The corrosponding row value for max col
%    
%     [yy2, ind] = min(y2); %find the min col in component n+1
%     xx2 = x2(ind); % The corrosponding row value for min col
%     
% %     Draw a line between the two points (xx1,yy1) and (xx2,yy2) and insert
% %     it in newboundaries4
%     if abs(yy1-yy2) < 50 
%         shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
%         newboundaries = step(shapeInserter, newboundaries, uint16([yy1 xx1 yy2 xx2]));
%     end
% %     figure, imshow(newboundaries4), title('After step shapeinserter');
%     
% end
% clear xx2 xx1 yy1 yy2 y1 y2 x1 x2;

figure, imshow(I,[]), title('Middle Connections')
%blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
%Use our diff1 as the AlphaData for the solid red image. 
set(displ, 'AlphaData', newboundaries)


%% Part 4, Logedges

zm_7_logedges;

% % clear 
% % close all
% % 
% % If not using rest of code, uncomment this part so this section can be run
% % separately!
% 
% % ptID = input('Enter image name you want to open: ','s'); %Request patient image name
% % ptID = strcat(ptID,'.tif'); 
% % I = imread(ptID); %open the image, keeping it in 16-bits
% % 
% % figure, imshow(I,[]) %to help decide if it should be cropped or not
% % title('Image 1: Original Image')
% 
% % detect edges using laplacian of gaussian method
% 
% figure;
% imshow(I,[]);
% title('Log Edges');
% magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %magenta has RGB value of 1 0 1
% hold on 
% h3 = imshow(magenta); 
% set(h3, 'AlphaData', BW_long);
% hold off
% 
% %overlapping final and log edges images
% 
% logfin = zeros(img_y,img_x);
% for al = 1:img_y
%     for bl = 1:img_x
%         if newboundaries(al,bl)==1
%             logfin(al,bl)=1;
%         end
%         if BW_long(al,bl)==1
%             logfin(al,bl)=1;
%         end
%     end
% end
% 
% figure;
% imshow(I,[]);
% title('Log and Final');
% green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
% hold on
% h4 = imshow(green);
% set(h4, 'AlphaData', logfin);
% hold off
% 
% % skeleton then connect lines
% connected = bwmorph(logfin,'skel',Inf);
% connected = bwmorph(connected,'bridge');
% connected = bwmorph(connected,'thicken');
% connected = bwareaopen(connected,40);
% %connected = bwmorph(connected,'clean');
% 
% figure;
% imshow(I,[]);
% title('Log and Final c,s');
% green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
% hold on
% h5 = imshow(green);
% set(h5, 'AlphaData', connected);
% hold off
% 
% %% Part 5, Final
% 
% connected=bwmorph(connected,'fill');
% connected1=bwmorph(connected,'bridge');
% connected2=bwmorph(connected1,'close');
% 
% 
% biggest = bwareafilt(connected2,1,'largest');
% 
% 
% 
% figure, imshow(I,[]), title('Thick Under Curve')
% %blue on top on figure
% blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
% hold on
% displ = imshow(blue);
% hold off
% set(displ, 'AlphaData', biggest)
% 
% bibi=bwmorph(biggest,'skel',Inf);
% figure, imshow(I,[]), title('Thin Under Curve')
% %blue on top on figure
% blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
% hold on
% displ = imshow(blue);
% hold off
% set(displ, 'AlphaData', bibi)
% lemmon=bibi;
% 
% %% Part 6, the Top Line
% 
% [r,c]=find(bibi == 1);
% 
% %makes a matrix with the columns in the first column and the rows in the
% %second
% matr(:,1)=c;
% matr(:,2)=r;
% 
% sortedbycol=sortrows(matr);
% rightside1=sortedbycol(1:round((1/2)*length(sortedbycol)),:);
% minright=sortrows(rightside1,2);
% leftside1=sortedbycol(round((1/2)*length(sortedbycol)):end,:);
% minleft=sortrows(leftside1,2);
% 
% 
% left=minleft(1,:);
% right=minright(1,:);
% 
% x1=right(1);
% y1=right(2);
% x2=left(1);
% y2=left(2);
% 
% %inserts the line across
% 
% shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
% check = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));
% 
% 
% shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
% connectedtop = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));
% 
% summed=0;
% [ii,kk]=size(bibi);
% for iii = 1:ii
%     for kkk = 1:kk
%         if check(iii,kkk)==1&& bibi(iii,kkk)==0
%             summed=summed+1;
%         end
%     end
% end
% 
% if summed<kk*(2/3)
%     disp('Single Connected Top Line Invalid')
% 
%     rightside=sortedbycol(1:round((1/3)*length(sortedbycol)),:);
%     minright=sortrows(rightside,2);
%     leftside=sortedbycol(round((2/3)*length(sortedbycol)):end,:);
%     minleft=sortrows(leftside,2);
%     midside=sortedbycol(round((1/3)*length(sortedbycol)):round((2/3)*length(sortedbycol)),:);
%     minmid=sortrows(midside,2);
% 
% 
%     left=minleft(1,:);
%     right=minright(1,:);
%     middlecol=minmid(1,:);
% 
%     x1=right(1);
%     y1=right(2);
%     x2=middlecol(1);
%     y2=middlecol(2);
%     x3=left(1);
%     y3=left(2);
%     
%     
%     shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
%     check = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));
% 
%     shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
%     fin = step(shapeInserter, check, uint16([x2 y2 x3 y3]));
% 
%     figure, imshow(I,[]), title('Connect Tops')
%     %blue on top on figure
%     blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
%     hold on 
%     displ = imshow(blue); 
%     hold off 
%     %Use our diff1 as the AlphaData for the solid red image. 
%     set(displ, 'AlphaData', fin)
% 
% 
% else
%     figure, imshow(I,[]), title('Connect Tops')
%     %blue on top on figure
%     blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
%     hold on 
%     displ = imshow(blue); 
%     hold off 
%     %Use our diff1 as the AlphaData for the solid red image. 
%     set(displ, 'AlphaData', connectedtop)
% end
% % 
% % 
% % BW2= imfill(connectedtop,'holes');
% % figure,imshow(BW2)
% % title('Filled in Chosen Space')
% % 
% % %% Part 7, masking out the rest of the space
% % finalimage=zeros(img_y,img_x);
% % darkestvalue=min(min(I));
% % for rt=1:img_x
% %     for nt=1:img_y
% %         if BW2(nt,rt)==1
% %             finalimage(nt,rt)=I(nt,rt);
% %         end
% %         if BW2(nt,rt)==0
% %             finalimage(nt,rt)=darkestvalue;
% %         end
% %     end
% % end
% % 
% % figure,imshow(finalimage,[]);
% % title('Only Breast Images')
