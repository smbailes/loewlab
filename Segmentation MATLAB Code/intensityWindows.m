 % split image in half
    mid_col = zeros(img_y,img_x);
    mid_col(:,img_x/2) = 1;

    %find coordinates of midline
    [midx midy] = find(mid_col==1); 
    midlinez(:,1) = midx;
    midlinez(:,2) = midy;

%     %find intersection between upper bound and midline
%     fprintf('Select Upper Bound\n');
%     figure, imshow(I, []), title('Bound Detection')
%     [Xup,Yup] = ginput(1);
%     ptofcompx = find(Yup==midlinez(:,1));
%     ptofcompy = midlinez(ptofcompx,2);
%     ptofcomp(:,1) = ptofcompx;
%     ptofcomp(:,2) = ptofcompy;
    
    %selects nipples 
    figure, imshow(I, []), title('Nipple Detection - Left')
    [Xleft,Yleft] = ginput(1);
    figure, imshow(I, []), title('Nipple Detection - Right')
    [Xright,Yright] = ginput(1);
    
    intwind = [];

    %find x and y distances of each pixel from ptofcomp and save into first 2 columns
    %use distances to find angle and save into 3rd column of matrix
    %find x1, y1, x2, and y2 for each pixel and save into columns 4-7
    [gettx getty] = find(circ_matrix==1);
    linepix1 = []; %creates matrix for coordinates on line 1
    linepix2 = []; %creates matrix for coordinates on line 2
    % len = 2; %length dimension for line - set at 2 FOR NOW (may need to change later) 
    for i = 1:length(gettx)
        if getty(i,1)>midy(1,1)
            intwind(i,1) = abs(Xleft - gettx(i,1));
            intwind(i,2) = abs(Yleft - getty(i,1));
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
        else
            intwind(i,1) = abs(Xright - gettx(i,1));
            intwind(i,2) = abs(Yright - getty(i,1));
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
    end 

    %find pixel intensities to corresponding coordinates 
    %should be 2 182x3 matrix (P10) 
    [ro co] = size(I);
    for i = 1:length(gettx)
        for j = 1:10
            xbel(j) = linepix1(i,j,1);
            if xbel(j)<1, xbel(j) = 1; end
            if xbel(j)>ro, xbel(j) = ro; end
            ybel(j) = linepix1(i,j,2);
            if ybel(j)<1, ybel(j) = 1; end
            intensities1(i,j) = I(round(xbel(j)),round(ybel(j))); 
            xab(j) = linepix2(i,j,1);
            if xab(j)<1, xab(j) = 1; end
            if xab(j)>ro, xab(j) = ro; end
            yab(j) = linepix2(i,j,2);
            if yab(j)<1, yab(j) = 1; end
            intensities2(i,j) = I(round(xab(j)),round(yab(j))); 
        end
        mean1 = mean(intensities1(i,1:j));
        mean2 = mean(intensities2(i,1:j));
        if abs(mean1-mean2)<10
            xrem = gettx(i);
            yrem = getty(i);
            circ_matrix(xrem,yrem) = 0; %removes pixels 
        end
    end 