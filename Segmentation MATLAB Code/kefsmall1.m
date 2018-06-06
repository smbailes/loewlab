%% Small Breasts Part 1, Circles
rmin = 75;
rmax = 125;
[centers, radii] = imfindcircles(edgecanny,[rmin rmax],'Sensitivity',0.97); 
%A higher 'Sensitivity' value sets the detection threshold lower and leads 
%to detecting more circles. 

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


    [e,f]=size(circ_matrix);
    lowers=zeros(e,f);
    for fir=1:f
        y=find(circ_matrix(:,fir)==1);
        if ~isempty(y)
        lowers(y(end),fir)=1;
        end
    end


    [dist1,dind1]=bwdist(lowers);
    gettem=find(dist1<2&dist1>0);
    if gettem==0
        gettem=1;
    end
    lowers(gettem)=1;


    figure, imshow(I,[]), title('Bottom of Circles')
    % magenta on top on figure
    magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %yellow has RGB value 1 1 0
    hold on 
    displ = imshow(magenta); 
    hold off 
    % Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', lowers) 

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


    [e,f]=size(circ_matrix);
    lowers=zeros(e,f);
    for fir=1:f
        y=find(circ_matrix(:,fir)==1);
        if ~isempty(y)
        lowers(y(end),fir)=1;
        end
    end


    [dist1,dind1]=bwdist(lowers);
    gettem=find(dist1<2&dist1>0);
    if gettem==0
        gettem=1;
    end
    lowers(gettem)=1;


    figure, imshow(I,[]), title('Bottom of Circles')
    % magenta on top on figure
    magenta = cat(3, ones(size(I)), zeros(size(I)), ones(size(I))); %yellow has RGB value 1 1 0
    hold on 
    displ = imshow(magenta); 
    hold off 
    % Use our diff1 as the AlphaData for the solid red image. 
    set(displ, 'AlphaData', lowers) 

    kefsmall2(Xup, Xlo, Yup, Ylo)
end

