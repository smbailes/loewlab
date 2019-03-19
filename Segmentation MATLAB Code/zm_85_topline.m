%% Top Line
% Places top line overlayed on final segmented image
% Last Updated by Zainab Mahmood on 2/5/18


%% Part 5, Final

connected=bwmorph(connected,'fill');
connected1=bwmorph(connected,'bridge');
connected2=bwmorph(connected1,'close');

[ro co] = size(I);

biggest = bwareafilt(connected2,1,'largest');
[bigarrayx bigarrayy] = find(biggest==1);
maxbig = max(bigarrayy); minbig = min(bigarrayy); 
count = 2;
while abs(maxbig-minbig) < ((3/4)*co)
    biggest = bwareafilt(connected2,count,'largest');
    [bigarrayx bigarrayy] = find(biggest==1);
    maxbig = max(bigarrayy); minbig = min(bigarrayy); 
    count = count+1;
end

% % % %     shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
% % % %     
% % % %     xxnew = []; yynew = [];
% % % %     [xxnew yynew] = find(biggest(:,1:(co/2)-1)==1); 
% % % %     yy1 = max(yynew);
% % % %     xx1 = xxnew(max(find(yynew == yy1))); 
% % % % 
% % % %     xxnew = []; yynew = [];
% % % %     [xxnew yynew] = find(biggest(:,(co/2):co)==1); 
% % % %     yy2 = min(yynew);
% % % %     xx2 = xxnew(max(find(yynew == yy2))); 
% % % %     yy2 = min(yynew)+(co/2);
% % % % 
% % % %     biggest = step(shapeInserter, biggest, uint16([yy1 xx1  yy2 xx2]));
    
% biggest(1:Yup, :) = 0;
% biggest(Ylo:ro, :) = 0;
% biggest(:,1:Xleft) = 0;
% biggest(:,Xright:co) = 0;


% % % % % % % % % % % [rconn cconn] = size(biggest);
% % % % % % % % % % % bnds = [colr coll];
% % % % % % % % % % % rowbnd = max(bnds);
% % % % % % % % % % % biggest(rowbnd:rconn, :) = 0; 
%biggest = connectDots(biggest,25);

figure, imshow(I,[]), title('Thick Under Curve')
%blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on
displ = imshow(blue);
hold off
set(displ, 'AlphaData', biggest)

bibi=bwmorph(biggest,'skel',Inf);
bibi(rlb:yb, :)=0;
figure, imshow(I,[]), title('Thin Under Curve')
%blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on
displ = imshow(blue);
hold off
set(displ, 'AlphaData', bibi)
lemmon=bibi;

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
%     figure(38), imshow(I,[]), title('Connect Tops')
%     %blue on top on figure
%     blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
%     hold on 
%     displ = imshow(blue); 
%     hold off 
%     %Use our diff1 as the AlphaData for the solid red image. 
%     set(displ, 'AlphaData', fin)
%     connectedtop=fin;
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

%% Manual Topline 

fprintf('Pick the top line of breast region. \n');
figure, imshow(I, [])


[rowtop coltop] = ginput(1);
coltop = round(coltop);
connected2(coltop, :) = 1; 
connectedtop = connected2;

figure, imshow(I,[]), title('Connect Tops')
%% Part 6, the Top Line

% % [r,c]=find(bibi == 1);

connectedtop = bibi; 

%makes a matrix with the columns in the first column and the rows in the
%second
% % matr(:,1)=c;
% % matr(:,2)=r;
% % 
% % sortedbycol=sortrows(matr);
% % rightside1=sortedbycol(1:round((1/2)*length(sortedbycol)),:);
% % minright=sortrows(rightside1,2);
% % leftside1=sortedbycol(round((1/2)*length(sortedbycol)):end,:);
% % minleft=sortrows(leftside1,2);
% % 
% % 
% % left=minleft(1,:);
% % right=minright(1,:);
% % 
% % x1=right(1);
% % y1=right(2);
% % x2=left(1);
% % y2=left(2);
% % 
% % %inserts the line across
% % 
% % shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
% % check = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));
% % 
% % 
% % shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
% % connectedtop = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));
% % 
% % summed=0;
% % [ii,kk]=size(bibi);
% % for iii = 1:ii
% %     for kkk = 1:kk
% %         if check(iii,kkk)==1&& bibi(iii,kkk)==0
% %             summed=summed+1;
% %         end
% %     end
% % end
% % 
% % if summed<kk*(2/3)
% %     disp('Single Connected Top Line Invalid')
% % 
% %     rightside=sortedbycol(1:round((1/3)*length(sortedbycol)),:);
% %     minright=sortrows(rightside,2);
% %     leftside=sortedbycol(round((2/3)*length(sortedbycol)):end,:);
% %     minleft=sortrows(leftside,2);
% %     midside=sortedbycol(round((1/3)*length(sortedbycol)):round((2/3)*length(sortedbycol)),:);
% %     minmid=sortrows(midside,2);
% % 
% % 
% %     left=minleft(1,:);
% %     right=minright(1,:);
% %     middlecol=minmid(1,:);
% % 
% %     x1=right(1);
% %     y1=right(2);
% %     x2=middlecol(1);
% %     y2=middlecol(2);
% %     x3=left(1);
% %     y3=left(2);
% %     
% %     
% %     shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
% %     check = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));
% % 
% %     shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
% %     fin = step(shapeInserter, check, uint16([x2 y2 x3 y3]));
% % 
% %     figure(38), imshow(I,[]), title('Connect Tops')
% %     %blue on top on figure
% %     blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
% %     hold on 
% %     displ = imshow(blue); 
% %     hold off 
% %     %Use our diff1 as the AlphaData for the solid red image. 
% %     set(displ, 'AlphaData', fin)
% %     connectedtop=fin;
% % 
% % else
    figure, imshow(I, []), title('Select Top Line')
    [Ytop,Xtop] = ginput(1);
    Xtop = round(Xtop);
    connectedtop(Xtop,:)=1;
    connectedtop(1:(Xtop-1),:)=0;
    
    figure, imshow(I,[]), title('Connect Tops')
    %blue on top on figure
    blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
    hold on 
    displ = imshow(blue); 
    hold off 
    %Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', connectedtop)
% % end
%% Finding Enclosed Area
% [rowc colc] = size(connectedtop);
% i = 1;
% r1 = 1; r2 = 1; c1 = 1; c2 = 1; 
% for a = 1:colc-1
%     locs = [];
%     for b = 1:rowc-1
%         if connectedtop(b,a)
%             if isempty(locs)
%                 r1 = b;
%                 c1 = a;
%             end
%             locs(i,1) = b;
%             locs(i,2) = a;
%             r2 = b;
%             c2 = a;
%             i = i + 1; 
%         end
%     end
%     %if r1~=0 && r2~=0 && r1~=r2
%         connectedtop(r1:r2, c1) = 1; 
%     %end
% end

%newCrop = imcrop(connectedtop);
%imshow(newCrop, [min(I) max(I)]);

% Invert the image to form black curves on a white background
connectedtop = ~connectedtop;
% Get rid of huge background that touches the border
connectedtop = imclearborder(connectedtop,4);
% Display the final image.
imshow(connectedtop, []);
title('Final Binary Image');

%% Cropping image
connectedtop=bwmorph(connectedtop,'fill');
connectedtop=bwmorph(connectedtop,'bridge');
connectedtop=bwmorph(connectedtop,'close');

figure, imshow(I,[]), title('Final Crop')
    %blue on top on figure
    blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
    hold on 
    displ = imshow(blue); 
    hold off 
    %Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', connectedtop)

imwrite(connectedtop, 'IRVT026-0015-CP.tif'); 

%% 

% % Create a binary image ("mask") from the ROI object.
% binaryImage = createMask(connectedtop);
% 
% %% 
% 
% % Get coordinates of the boundary of the freehand drawn region.
% structBoundaries = bwboundaries(binaryImage);
% xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
% x = xy(:, 2);   % Columns.
% y = xy(:, 1);   % Rows.
% % Mask the image outside the mask, and display it.
% 
% % Will keep only the part of the image that's inside the mask, zero outside mask.
% blackMaskedImage = firstImage;
% blackMaskedImage(~binaryImage) = 0;
% 
% % Now crop the image.
% leftColumn = min(x);
% rightColumn = max(x);
% topLine = min(y);
% bottomLine = max(y);
% width = rightColumn - leftColumn + 1;
% height = bottomLine - topLine + 1;    
% 
% newCrop = imcrop(blackMaskedImage, [leftColumn, topLine, width, height]);
% close;
% % 
% 
% BW2= imfill(connectedtop,'holes');
% figure,imshow(BW2)
% title('Filled in Chosen Space')
% 
% %% Part 7, masking out the rest of the space
% finalimage=zeros(img_y,img_x);
% darkestvalue=min(min(I));
% for rt=1:img_x
%     for nt=1:img_y
%         if BW2(nt,rt)==1
%             finalimage(nt,rt)=I(nt,rt);
%         end
%         if BW2(nt,rt)==0
%             finalimage(nt,rt)=darkestvalue;
%         end
%     end
% end
% 
% figure,imshow(finalimage,[]);
% title('Only Breast Images')
