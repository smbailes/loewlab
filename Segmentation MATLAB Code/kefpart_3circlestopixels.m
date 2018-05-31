circlines=circles.Children;

[img_y, img_x] = size(I);

      %create matrix of zeros as big as the image
        xc1 = circlines(1).XData;
        yc1 = circlines(1).YData;

xc1(isnan(xc1)) = [];
yc1(isnan(yc1)) = [];

if length(centers)==2
    firsx=xc1(1:(length(xc1)/2));
    firsy=yc1(1:(length(yc1)/2));

    lasx=xc1((length(xc1)/2)+1:end);
    lasy=yc1((length(yc1)/2)+1:end);

    circ_matrix1=zeros(img_y,img_x);
    for d = 1:length(firsx)
        firsx(d) = round(firsx(d));
        firsy(d) = round(firsy(d));
    if firsy(d)<0.5
        firsy(d)=1;
    end
    if firsx(d)<0.5
        firsx(d)=1;
    end
        circ_matrix1(firsy(d),firsx(d)) = 1;      %fill in 1's wherever there is a point in the circle
    end

    circ_matrix2=zeros(img_y,img_x);
    for d = 1:length(lasx)
        lasx(d) = round(lasx(d));
        lasy(d) = round(lasy(d));
    if lasy(d)<0.5
        lasy(d)=1;
    end
    if lasx(d)<0.5
        lasx(d)=1;
    end
        circ_matrix2(lasy(d),lasx(d)) = 1;      %fill in 1's wherever there is a point in the circle
    end
   

    figure, imshow(I,[]), title('Right Circle Pixels')
    % magenta on top on figure
    magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %yellow has RGB value 1 1 0
    hold on 
    displ = imshow(magenta); 
    hold off 
    % Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', circ_matrix1) 

    figure, imshow(I,[]), title('Left Circle Pixels')
    % green on top on figure
    green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I))); %yellow has RGB value 1 1 0
    hold on 
    displ = imshow(green); 
    hold off 
    % Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', circ_matrix2)   

end

circ_matrix=zeros(img_y,img_x);
for d = 1:length(xc1)
    xc1(d) = round(xc1(d));
    yc1(d) = round(yc1(d));
    if yc1(d)<0.5
        yc1(d)=1;
    end
    if xc1(d)<0.5
        xc1(d)=1;
    end
    circ_matrix(yc1(d),xc1(d)) = 1;      %fill in 1's wherever there is a point in the circle
end
    
    
figure, imshow(I,[]), title('Both Circle Pixels')
% yellow on top on figure
yellow = cat(3, ones(size(I)), ones(size(I)), zeros(size(I))); %yellow has RGB value 1 1 0
hold on 
displ = imshow(yellow); 
hold off 
% Use our diff1 as the AlphaData for the solid red image. 
set(displ, 'AlphaData', circ_matrix) 


%finding bottom half of (both) circles:
[e,f]=size(circ_matrix);
circok=zeros(e,f);
for fir=1:f
    y=find(circ_matrix(:,fir)==1);
    if ~isempty(y)
    circok(y(end),fir)=1;
    end
end

okay=bwmorph(circok,'bridge');
okay2=bwmorph(okay,'thicken');
okay3=bwmorph(okay2,'bridge');
circs=bwmorph(okay3,'thicken');

figure, imshow(I,[]), title('Lower Half Overlapped Pixels')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
% Use our diff1 as the AlphaData for the solid red image. 
set(displ, 'AlphaData', circs)


in = input('Is the breast small or large? Enter s/l: ','s');
if in == 's'
    kefpart_5SB_Points
elseif in == 'l'
    kefpart_4Points
end