%% Top Line
% Places top line overlayed on final segmented image
% Last Updated by Zainab Mahmood on 2/5/18


%% Part 5, Final

connected=bwmorph(connected,'fill');
connected1=bwmorph(connected,'bridge');
connected2=bwmorph(connected1,'close');


biggest = bwareafilt(connected2,1,'largest');



figure, imshow(I,[]), title('Thick Under Curve')
%blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on
displ = imshow(blue);
hold off
set(displ, 'AlphaData', biggest)

bibi=bwmorph(biggest,'skel',Inf);
figure, imshow(I,[]), title('Thin Under Curve')
%blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on
displ = imshow(blue);
hold off
set(displ, 'AlphaData', bibi)
lemmon=bibi;

%% Part 6, the Top Line

[r,c]=find(bibi == 1);

%makes a matrix with the columns in the first column and the rows in the
%second
matr(:,1)=c;
matr(:,2)=r;

sortedbycol=sortrows(matr);
rightside1=sortedbycol(1:round((1/2)*length(sortedbycol)),:);
minright=sortrows(rightside1,2);
leftside1=sortedbycol(round((1/2)*length(sortedbycol)):end,:);
minleft=sortrows(leftside1,2);


left=minleft(1,:);
right=minright(1,:);

x1=right(1);
y1=right(2);
x2=left(1);
y2=left(2);

%inserts the line across

shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
check = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));


shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
connectedtop = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));

summed=0;
[ii,kk]=size(bibi);
for iii = 1:ii
    for kkk = 1:kk
        if check(iii,kkk)==1&& bibi(iii,kkk)==0
            summed=summed+1;
        end
    end
end

if summed<kk*(2/3)
    disp('Single Connected Top Line Invalid')

    rightside=sortedbycol(1:round((1/3)*length(sortedbycol)),:);
    minright=sortrows(rightside,2);
    leftside=sortedbycol(round((2/3)*length(sortedbycol)):end,:);
    minleft=sortrows(leftside,2);
    midside=sortedbycol(round((1/3)*length(sortedbycol)):round((2/3)*length(sortedbycol)),:);
    minmid=sortrows(midside,2);


    left=minleft(1,:);
    right=minright(1,:);
    middlecol=minmid(1,:);

    x1=right(1);
    y1=right(2);
    x2=middlecol(1);
    y2=middlecol(2);
    x3=left(1);
    y3=left(2);
    
    
    shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
    check = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));

    shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
    fin = step(shapeInserter, check, uint16([x2 y2 x3 y3]));

    figure, imshow(I,[]), title('Connect Tops')
    %blue on top on figure
    blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
    hold on 
    displ = imshow(blue); 
    hold off 
    %Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', fin)


else
%     figure, imshow(I,[]), title('Connect Tops')
%     %blue on top on figure
%     blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
%     hold on 
%     displ = imshow(blue); 
%     hold off 
%     %Use our diff1 as the AlphaData for the solid red image. 
%     set(displ, 'AlphaData', connectedtop)
end
% 
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
