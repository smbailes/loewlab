%% Small Breasts Part 1, Circles
rmin = 75;%75
rmax = 250; %125
[centers, radii] = imfindcircles(edgecanny3,[rmin rmax],'Sensitivity',0.97); 
%A higher 'Sensitivity' value sets the detection threshold lower and leads 
%to detecting more circles. 
%% 


[cr,cc]=size(centers);

if cr<2
    imshow(I,[]);
    hold on
    circle1 = viscircles(centers(1,:),radii(1,:),'LineWidth',1);
    
    circlines=circle1.Children;
    xc1 = circlines(1).XData;
    yc1 = circlines(1).YData;
    xc1(isnan(xc1)) = [];
    yc1(isnan(yc1)) = [];
    xc1=round(xc1);
    yc1=round(yc1);

    [img_y,img_x]=size(I);
    circ_matrix=zeros(img_y,img_x);
    for d=1:length(xc1)
        circ_matrix(yc1(d),xc1(d)) = 1;
    end

    figure, imshow(I,[]), title('Circles')
    % magenta on top on figure
    magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %yellow has RGB value 1 1 0
    hold on 
    displ = imshow(magenta); 
    hold off 
    % Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', circ_matrix) 
    
    % split image in half
    mid_col = zeros(img_y,img_x);
    mid_col(:,img_x/2) = 1;

    %find coordinates of midline
    [midx midy] = find(mid_col==1); 
    midlinez(:,1) = midx;
    midlinez(:,2) = midy;

    %find intersection between upper bound and midline
    ptofcompx = find(Yup==midlinez(:,1));
    ptofcompy = midlinez(ptofcompx,2);
    ptofcomp(:,1) = ptofcompx;
    ptofcomp(:,2) = ptofcompy;

    intwind = [];

    %find x and y distances of each pixel from ptofcomp and save into first 2 columns
    %use distances to find angle and save into 3rd column of matrix
    %find x1, y1, x2, and y2 for each pixel and save into columns 4-7
    [gettx getty] = find(circ_matrix==1);
    linepix1 = []; %creates matrix for coordinates on line 1
    linepix2 = []; %creates matrix for coordinates on line 2
    % len = 2; %length dimension for line - set at 2 FOR NOW (may need to change later) 
    for i = 1:length(gettx)
        intwind(i,1) = abs(ptofcompx - gettx(i,1));
        intwind(i,2) = abs(ptofcompy - getty(i,1));
        intwind(i,3) = atan(intwind(i,2)/intwind(i,1));
        for j = 1:10 %line dimension set as 3 pixels
            xdiff(i,j) = j*cos(intwind(i,3));
            ydiff(i,j) = j*sin(intwind(i,3));
        end
        for k = 1:10
            linepix1(i,k,1) = intwind(i,1) - xdiff(i,k);
            linepix1(i,k,2) = intwind(i,1) - xdiff(i,k);
            linepix2(i,k,1) = intwind(i,1) + xdiff(i,k);
            linepix2(i,k,2) = intwind(i,1) + xdiff(i,k);
        end
    end 

    %find pixel intensities to corresponding coordinates 
    %should be 2 182x3 matrix (P10) 
    for i = 1:length(gettx)
        for j = 1:10
            xbel(j) = linepix1(i,j,1);
            ybel(j) = linepix1(i,j,2);
            intensities1(i,j) = I(round(xbel(j)),round(ybel(j))); 
            xab(j) = linepix2(i,j,1);
            yab(j) = linepix2(i,j,2);
            intensities2(i,j) = I(round(xab(j)),round(yab(j))); 
        end
        mean1 = mean(intensities1(i,1:j));
        mean2 = mean(intensities2(i,1:j));
        if abs(mean1-mean2)<100
            xrem = gettx(i);
            yrem = getty(i);
            circ_matrix(xrem,yrem) = 0; %removes pixels 
        end
    end 


%     [e,f]=size(circ_matrix);
%     lowers=zeros(e,f);
%     for fir=1:f
%         y=find(circ_matrix(:,fir)==1);
%         if ~isempty(y)
%         lowers(y(end),fir)=1;
%         end
%     end
% 
% 
%     [dist1,dind1]=bwdist(lowers);
%     gettem=find(dist1<2&dist1>0);
%     if gettem==0
%         gettem=1;
%     end
%     lowers(gettem)=1;


    figure, imshow(I,[]), title('Bottom of Circles')
    % magenta on top on figure
    magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %yellow has RGB value 1 1 0
    hold on 
    displ = imshow(magenta); 
    hold off 
    % Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', circ_matrix) 

    kefsmall2

else 
	if cr>2   % keeps it to only 2 circles. (hopefully the correct ones)
        centers(3:end,:)=[];
        radii(3:end,:)=[];
	end
	imshow(I,[]);
    hold on
    circle1 = viscircles(centers(1,:),radii(1,:),'LineWidth',1);
    circle2 = viscircles(centers(2,:),radii(2,:),'LineWidth',1);

    circlines=circle1.Children;
    xc1 = circlines(1).XData;
    yc1 = circlines(1).YData;
    xc1(isnan(xc1)) = [];
    yc1(isnan(yc1)) = [];
    xc1=ceil(xc1);
    yc1=ceil(yc1);
    for r=1:length(xc1)
        if xc1(r)==0
            xc1(r)=1;
        end
        if yc1(r)==0
            yc1(r)=1;
        end
    end

    circlines2=circle2.Children;
    xc2 = circlines2(1).XData;
    yc2 = circlines2(1).YData;
    xc2(isnan(xc2)) = [];
    yc2(isnan(yc2)) = [];
    xc2=ceil(xc2);
    yc2=ceil(yc2);
    for r=1:length(xc2)
        if xc2(r)==0
            xc2(r)=1;
        end
        if yc2(r)==0
            yc2(r)=1;
        end
    end
    
    [img_y,img_x]=size(I);
    circ_matrix=zeros(img_y,img_x);
    for d=1:length(xc1)
        circ_matrix(yc1(d),xc1(d)) = 1;
        circ_matrix(yc2(d),xc2(d)) = 1;
    end

    figure, imshow(I,[]), title('Circles')
    % magenta on top on figure
    magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %yellow has RGB value 1 1 0
    hold on 
    displ = imshow(magenta); 
    hold off 
    % Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', circ_matrix) 

    % split image in half
    mid_col = zeros(img_y,img_x);
    mid_col(:,img_x/2) = 1;

    %find coordinates of midline
    [midx midy] = find(mid_col==1); 
    midlinez(:,1) = midx;
    midlinez(:,2) = midy;

    %find intersection between upper bound and midline
    fprintf('Select Upper Bound\n');
    figure, imshow(I, []), title('Bound Detection')
    [Xup,Yup] = ginput(1);
    ptofcompx = find(Yup==midlinez(:,1));
    ptofcompy = midlinez(ptofcompx,2);
    ptofcomp(:,1) = ptofcompx;
    ptofcomp(:,2) = ptofcompy;

    intwind = [];

    %find x and y distances of each pixel from ptofcomp and save into first 2 columns
    %use distances to find angle and save into 3rd column of matrix
    %find x1, y1, x2, and y2 for each pixel and save into columns 4-7
    [gettx getty] = find(circ_matrix==1);
    linepix1 = []; %creates matrix for coordinates on line 1
    linepix2 = []; %creates matrix for coordinates on line 2
    % len = 2; %length dimension for line - set at 2 FOR NOW (may need to change later) 
    for i = 1:length(gettx)
        intwind(i,1) = abs(ptofcompx - gettx(i,1));
        intwind(i,2) = abs(ptofcompy - getty(i,1));
        intwind(i,3) = atan(intwind(i,2)/intwind(i,1));
        for j = 1:10 %line dimension set as 3 pixels
            xdiff(i,j) = j*cos(intwind(i,3));
            ydiff(i,j) = j*sin(intwind(i,3));
        end
        for k = 1:10
            linepix1(i,k,1) = intwind(i,1) - xdiff(i,k);
            linepix1(i,k,2) = intwind(i,1) - xdiff(i,k);
            linepix2(i,k,1) = intwind(i,1) + xdiff(i,k);
            linepix2(i,k,2) = intwind(i,1) + xdiff(i,k);
        end
    end 

    %find pixel intensities to corresponding coordinates 
    %should be 2 182x3 matrix (P10) 
    for i = 1:length(gettx)
        for j = 1:10
            xbel(j) = linepix1(i,j,1);
            ybel(j) = linepix1(i,j,2);
            intensities1(i,j) = I(round(xbel(j)),round(ybel(j))); 
            xab(j) = linepix2(i,j,1);
            yab(j) = linepix2(i,j,2);
            intensities2(i,j) = I(round(xab(j)),round(yab(j))); 
        end
        mean1 = mean(intensities1(i,1:j));
        mean2 = mean(intensities2(i,1:j));
        if abs(mean1-mean2)<100
            xrem = gettx(i);
            yrem = getty(i);
            circ_matrix(xrem,yrem) = 0; %removes pixels 
        end
    end 

%     [e,f]=size(circ_matrix);
%     lowers=zeros(e,f);
%     for fir=1:f
%         y=find(circ_matrix(:,fir)==1);
%         if ~isempty(y)
%         lowers(y(end),fir)=1;
%         end
%     end
% 
% 
%     [dist1,dind1]=bwdist(lowers);
%     gettem=find(dist1<2&dist1>0);
%     if gettem==0
%         gettem=1;
%     end
%     lowers(gettem)=1;


    figure, imshow(I,[]), title('Bottom of Circles')
    % magenta on top on figure
    magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %yellow has RGB value 1 1 0
    hold on 
    displ = imshow(magenta); 
    hold off 
    % Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', circ_matrix) 

    kefsmall2
end

