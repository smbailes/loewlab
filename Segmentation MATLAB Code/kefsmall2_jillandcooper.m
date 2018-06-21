%%Small Points
close all;
smallpoints=zeros(img_y,img_x);

for aa = 1:img_y
    for bb = 1:img_x
        if edgecanny(aa,bb)==1
            smallpoints(aa,bb)=smallpoints(aa,bb)+2;
        end
        if lowers(aa,bb)==1
            smallpoints(aa,bb)=smallpoints(aa,bb)+3;
        end
        if newI(aa,bb)~=0
            smallpoints(aa,bb)=smallpoints(aa,bb)+1;
        end
%         if ellipses(aa,bb)==1
%             smallpoints(aa,bb)=smallpoints(aa,bb)+1;
%         end
    end
end


figure,imshow(smallpoints,[])
title('Combined Point Systems(S)')

smalloverlayedpoints = zeros(img_y,img_x);
for cc=1:img_y
    for dd=1:img_x
        if smallpoints(cc,dd)>1
            smalloverlayedpoints(cc,dd)=1;
        end
    end
end


[distover,dindover]=bwdist(smalloverlayedpoints);
righto=find(distover<2&distover>0);
if righto==0
    righto=1;
end
smalloverlayedpoints(righto)=1;


figure, imshow(I,[]), title('Small Pre-Clean')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', smalloverlayedpoints)

%% Part 2, Small Clean
sm=bwmorph(smalloverlayedpoints,'close');
sk=bwareaopen(sm,75);
sb=bwmorph(sk,'endpoints');
sl=bwmorph(sb,'clean');

[dist,dind]=bwdist(sl);


for a=1:length(dist) 
if dist(a)==1
     sl(dind(a))=1;
end
end


strell = strel('disk',2); %Create a Morphological structuring element, you change the shape used and diameter
sl= imclose(sl,strell); 


C2 = bwconncomp(sl); %Find connected components

long(2)=0; %preallocation
for f=1:C2.NumObjects
    long(f)=length(C2.PixelIdxList{f});
end
[num,ind]=max(long);
beforeconnect=zeros(img_y,img_x);
beforeconnect(C2.PixelIdxList{ind})=1;

C2.PixelIdxList{ind}=0;

seclong(2)=0;
for f=1:C2.NumObjects
    seclong(f)=length(C2.PixelIdxList{f});
end
[num2,ind2]=max(seclong);
beforeconnect(C2.PixelIdxList{ind2})=1;



figure, imshow(I,[]), title('removed')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on
displ = imshow(blue);
hold off
set(displ, 'AlphaData', beforeconnect)

[rowsbc,colsbc]=size(beforeconnect);
%% 
% Ellipses clean 

    %RIGHT SIDE
% override some default parameters
paramsr2.minMajorAxis = 150;
paramsr2.maxMajorAxis = 300;
paramsr2.numBest = 12; %draws 12 ellipses
paramsr2.rotation = 45; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
paramsr2.rotationSpan = 35;
%paramsr2.randomize = 0; %randomization component that may reduce changing of
%ellipses

% note that the edge (or gradient) image is used
bestFits2r = ellipseDetection(beforeconnect(:,1:round(colsbc/2)), paramsr2);
fprintf('Output %d best fits.\n', size(bestFits2r,1));


% ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 

% takes the information that was found of the ellipses and draws them;also
% keeping the information for each ellipse in a cell in qr(and later ql for those):
secondqr{1,length(bestFits2r)}=0;
for n=1:length(bestFits2r)
    secondqr{n} = ellipse(bestFits2r(n,3),bestFits2r(n,4),bestFits2r(n,5)*pi/180,bestFits2r(n,1),bestFits2r(n,2),'k');
end



%overriding parameters:
paramslr.minMajorAxis = 150;
paramslr.maxMajorAxis = 300;
paramslr.numBest = 12; %draws 12 ellipses
paramslr.rotation = 135; %If rotationSpan is in (0,90), only angles within [rotation-rotationSpan,rotation+rotationSpan] are accepted.
paramslr.rotationSpan = 43;
%paramslr.randomize = 0; %randomization component that may reduce changing of
%ellipses


%LEFT SIDE
bestFitslr = ellipseDetection(beforeconnect, paramslr);
fprintf('Output %d best fits.\n', size(bestFitslr,1));

%ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
secondql{1,length(bestFitslr)}=0;
for n=1:length(bestFitslr)
    secondql{n} = ellipse(bestFitslr(n,3),bestFitslr(n,4),bestFitslr(n,5)*pi/180,bestFitslr(n,1),bestFitslr(n,2),'k');
end


[img_y, img_x] = size(I);
laterellipses=zeros(img_y,img_x);

% draw in LEFT breast pixel by pixel
for a = 1:length(secondql)                %a,b,n,m are just used as counters in the for loops - delete at end of section
    e1 = secondql{a};
    for b = 1:length(e1.XData)              %get x and y data from cell array of ellipses and round so we can use them as indices
        xe1(a,b) = e1.XData(b);
        xe1(a,b) = round(xe1(a,b));
        ye1(a,b) = e1.YData(b);
        ye1(a,b) = round(ye1(a,b));
        if xe1(a,b)<0.5
            xe1(a,b)=1;
        end
        if ye1(a,b)<0.5
            ye1(a,b)=1;
        end
    end
    for d = 1:length(xe1)
        laterellipses(ye1(a,d),xe1(a,d)) = 1;      %fill in 1's wherever there is a point in the ellipse
    end
    
end

% draw in RIGHT breast pixel by pixel
for a = 1:length(secondqr)                %a,b,n,m are just used as counters in the for loops - delete at end of section
    e2 = secondqr{a};
    for b = 1:length(e2.XData)
        xe2(a,b) = e2.XData(b);
        xe2(a,b) = round(xe2(a,b));
        ye2(a,b) = e2.YData(b);
        ye2(a,b) = round(ye2(a,b));   
        if xe2(a,b)<0.5
            xe2(a,b)=1;
        end
        if ye2(a,b)<0.5
            ye2(a,b)=1;
        end
    end
    for d = 1:length(xe2)
        laterellipses(ye2(a,d),xe2(a,d)) = 1;      %fill in 1's wherever there is a point in the ellipse
    end
end

[checky,checkx]=size(laterellipses);
if checkx > img_x
    laterellipses=laterellipses(1:img_y,1:img_x);
end
if checky > img_y
    laterellipses=laterellipses(1:img_y,1:img_x);
end
    
[e,f]=size(laterellipses);

[elliprows,ellipcols]=find(laterellipses>0);
maxx=max(elliprows);
minn=min(elliprows);
diff=maxx-minn;
topthird=round((diff*2)/5);

laterellipses(1:topthird,:)=0;

% 
% for fir=1:f
%     y=find(laterellipses(:,fir)==1);
%     if ~isempty(y)
%     bott(y(end),fir)=1;
%     end
% end

figure, imshow(I,[]), title('Lower Half Second Ellipses')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on
displ = imshow(blue);
hold off
set(displ, 'AlphaData', laterellipses)
    
%% Part 3, Small Second set of points
secondpoints=zeros(img_y,img_x);

for aa = 1:img_y
    for bb = 1:img_x

        if lowers(aa,bb)==1
            secondpoints(aa,bb)=secondpoints(aa,bb)+3;
        end
        if newI(aa,bb)~=0
            secondpoints(aa,bb)=secondpoints(aa,bb)+1;
        end 
        if laterellipses(aa,bb)==2 %WAS 1!!!!
            secondpoints(aa,bb)=secondpoints(aa,bb)+2;
        end
    end
end

figure,imshow(smallpoints,[])
title('Second Points')

secondoverlap = zeros(img_y,img_x);
for cc=1:img_y
    for dd=1:img_x
        if secondpoints(cc,dd)>1
            secondoverlap(cc,dd)=1;
        end
    end
end


[distover,dindover]=bwdist(secondoverlap);
righto=find(distover<2&distover>0);
if righto==0
    righto=1;
end
secondoverlap(righto)=1;


figure, imshow(I,[]), title('Second Cleaning')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', secondoverlap)


secc=bwmorph(secondoverlap,'close');
figure, imshow(I,[]), title('Secc')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', secc)


skels=bwmorph(secc,'skel',Inf);
figure, imshow(I,[]), title('skels')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on
displ = imshow(blue);
hold off
set(displ, 'AlphaData', skels)


[rl,cl]=size(skels);

[xs,ys]=find(skels>0); %(1:round(rl/2),round(cl/3):round((cl*2)/3))



vend=zeros(rl,cl);


for s=1:length(xs)
    vend(xs(s),ys(s))=1;
end

vend(:,1:round((1/3)*cl))=0;
vend(:,round((2/3)*cl):end)=0;



figure, imshow(I,[]), title('vend')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on
displ = imshow(blue);
hold off
set(displ, 'AlphaData', vend)

gett=vend;


%% Part 4, 

CC = bwconncomp(gett);
newboundaries = gett;

for n = 1:CC.NumObjects - 1 %the number of lines in the middle region of the patient
%     Store all row and col values of component n and the component after in
%     x1,y1, x2, y2
    [x1, y1] = ind2sub(size(newboundaries),CC.PixelIdxList{n});
    [x2, y2] = ind2sub(size(newboundaries),CC.PixelIdxList{n+1});
    
    [yy1, ind] = min(x1); %find the max col in component n
    xx1 = y1(ind); % The corrosponding row value for max col
    
    [yy2, ind] = min(x2); %find the min col in component n+1
    xx2 = y2(ind); % The corrosponding row value for min col
    
    %Draw a line between the two points (xx1,yy1) and (xx2,yy2) and insert
    %it in newboundaries4
    shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
    newboundaries = step(shapeInserter, newboundaries, uint16([xx1 yy1 xx2 yy2]));
    %figure, imshow(newboundaries), title('After step shapeinserter');
    
end
clear xx2 xx1 yy1 yy2 y1 y2 x1 x2;

figure, imshow(I,[]), title('Middle Connections')
%blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
%Use our diff1 as the AlphaData for the solid red image. 
set(displ, 'AlphaData', newboundaries)



%% Part 5, Add again

finalbeforelog=zeros(img_y,img_x);
for go=1:img_y
    for rtw=1:img_x
        if newboundaries(go,rtw)==1
            finalbeforelog(go,rtw)=1;
        end
        if skels(go,rtw)==1
            finalbeforelog(go,rtw)=1;
        end
    end
end
newboundaries=finalbeforelog;


%% Part 6, zm_7_logedges_og

BW = edge(I,'log');

BW_long = bwareaopen(BW,20);
BW_long = bwmorph(BW_long,'thicken');

figure;
imshow(I,[]);
title('Log Edges');
magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %magenta has RGB value of 1 0 1
hold on 
h3 = imshow(magenta); 
set(h3, 'AlphaData', BW_long);
hold off

%overlapping final and log edges images

logfin = zeros(img_y,img_x);
for al = 1:img_y
    for bl = 1:img_x
        if newboundaries(al,bl)==1
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

% skeleton then connect lines
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

%% Part 7

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

%% Part 8, Finding top connector

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



shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White','LineWidth',1);
connectedtop = step(shapeInserter, bibi, uint16([x1 y1 x2 y2]));


%figure, imshow(I,[]), title('Connect Tops')
%blue on top on figure
%blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
%hold on 
%displ = imshow(blue); 
%hold off 
%Use our diff1 as the AlphaData for the solid red image. 
%set(displ, 'AlphaData', connectedtop)


%JILLIAN
%find the left most endpoint of breast outline
% (x1,y1)
%[x1,y1]=find()
%find the right most endpoint of the breast outline
% (x2,y2)
%{x2,y2]=find()
%create a line connecting the endpoints
%hold on
%line([x1 x2], [y1 y2])
figure, imshow(I,[]), title('Connect Tops')
hold on
BW2 = bwmorph(connectedtop,'skel',inf);
endPoints = bwmorph(BW2, 'endpoints');
%imshow(BW2);hold on;
imshowpair(I,BW2,'diff');hold on;    %imshowpair with 'diff' allows for the image to be unaltered
[rows cols] = find(endPoints);
[lowestX indexOfLowestX] = min(rows);
[highestX indexOfHighestX] = max(cols);
sonY = [rows(indexOfLowestX) rows(indexOfHighestX)];
sonX = [lowestX highestX];
plot(sonX, sonY, 'go');
hold on
%line([rows(indexOfLowestX) rows(indexOfHighestX)], [lowestx highestx])
horiz=line([0:highestX], [rows(indexOfLowestX)*ones(size(0:highestX))])
plot(0:highestX, rows(indexOfLowestX)*ones(size(0:highestX))) %highestx is x value of red, rows(indexOfHighestX) y coordinate
hold on
vert=line([highestX*ones(size(0:rows(indexOfHighestX)))],[0:rows(indexOfHighestX)])
plot(highestX*ones(size(0:rows(indexOfHighestX))),0:rows(indexOfHighestX))
hold on

xval1=[0:highestX];
xval2=[highestX*ones(size(0:rows(indexOfHighestX)))];
yval1=[rows(indexOfLowestX)*ones(size(0:highestX))];
yval2=[0:rows(indexOfHighestX)];

xvalintersect=find(xval1==xval2(1));
disp(xvalintersect)
yvalintersect=find(yval2==yval1(1));
disp(yvalintersect)

plot(xvalintersect,yvalintersect,'go');
hold off

n=0;
nn=0;
intersectpointx=1;
intersectpointy=1;

while intersectpointx~=0    
  while intersectpointy~=0
      ytest1=yval1(n+1);
      ytest2=yval2(n+1);
      intersectpointy=ytest2-ytest1;
      n=n+1;
  end
    xtest1=xval1(nn+1);
    xtest2=xval2(1);
    intersectpointx=xtest1-xtest2;
    nn=nn+1;
end



%ipoint=intersect(horiz,vert)
%plot(ipoint)
%plot(0:rows(indexOfHighestX), lowestX*ones(size(0:rows(indexOfHighestX))))
%line([0 sonX], [0 inf], 'rs')


% BW2= imfill(connectedtop,'holes');
% figure,imshow(BW2)
% title('Filled in Chosen Space')
% 
% %% Part 9, masking out the rest of the space
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
