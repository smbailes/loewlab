%% Clear Everything
clc;
close all; clear all;

%% Input and Show Image

snakeID = input('Enter image name you want to open: ','s'); %Request patient image name
snakeID = strcat(snakeID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Snakes Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
I = imread([dir '/' snakeID]); 

figure, imshow(I,[]) %to help decide if it should be cropped or not
title('Snakes Image')

%% Input Points and Setup Values

[yp,xp] = getpts;
Cur_We = 1 ;
Con_We = 1 ;
Gra_We = 1 ;
cur_Tresh = 0.8 ;
%con_Tresh = handles.metricdata.Con_Th ;
ws = 3 ;
max_iter = 30 ;
Stop_Cr = 0.0001 ;

Cur_We = Cur_We/(Cur_We + Con_We + Gra_We) ;
Con_We = Con_We/(Cur_We + Con_We + Gra_We) ;
Gra_We = Gra_We/(Cur_We + Con_We + Gra_We) ;

%% Run Code

gx = fspecial('sobel')/8 ;
gy = gx';
% cur_Tresh = 0.8 ;
% con_Tresh = 0.9 ; 
%I = imread('041_2_3.bmp') ;
max_I = max(I(:));
rn = size(I,1);
cn = size(I,2);
It = double(I);
Itt = It;

xs = length(xp);
ne1(1) = xs;
ne2(1) = 2;
for i = 2:xs-1
    ne1(i) = i-1;
    ne2(i) = i+1;
end
ne2(xs) = 1;
ne1(xs) = xs-1;    

gradx = imfilter(It,gx,'same');
grady = imfilter(It,gy,'same');
Igrad = (gradx.^2 + grady.^2).^(1/2);
% ws = 4 ;
Old_Etot = 100;
New_Etot = 1;
%for iter = 1:60
iter = 0;
Tuv = ones(1,xs);
Stop_Force = 0;
while(((abs(New_Etot - Old_Etot) > Stop_Cr) | (iter < 2)) & (iter < max_iter) & (Stop_Force == 0))
    iter = iter + 1;
    Old_Etot = New_Etot;
    Itt = It;
    lix = [];
    liy = [];
    for i = 1:xs
        r1 = max([1 xp(i)-3]):min([rn xp(i)+3]) ;
        r2 = max([1 yp(i)-3]):min([cn yp(i)+3]) ;
        zIt = Itt(fix(r1),fix(r2)) ;
        mean_Itt = mean(zIt(:)) ;
        max_Itt = max(zIt(:)) ;
        min_Itt = min(zIt(:)) ;
        lix = [lix [yp(i) ; yp(ne1(i))]] ;
        liy = [liy [xp(i) ; xp(ne1(i))]] ;
        if(~isempty(r1) & ~isempty(r2))
    %         It(fix(r2),fix(r1)) = ones(length(r2),length(r1))*255 ;
%            Itt(fix(r1),fix(r2)) = ones(length(r1),length(r2))*double((255 - (mean_Itt))) ;
            Itt(fix(r1),fix(r2)) = 255 - Itt(fix(r1),fix(r2)) ;
        else
            Tuv(i) = 0 ;
        end
    end
    pause(0.02),
    
    imshow(Itt,[])    
    line(lix,liy,'Color',[1 0 0],'linewidth',1) ;
    
    for i =1:xs
        ac = [xp(i) yp(i)];
        n1 = ne1(i) ;
        d(i) = ((xp(n1) - ac(1))^2 + (yp(n1) - ac(2))^2)^(1/2) ;
    end
    d_mean = mean(d) ;
    Econ_Te = cell(1,xs) ; 
    Ecur_Te = cell(1,xs) ; 
    Egra_Te = cell(1,xs) ; 
    c_con = zeros(1,xs) ;
    c_cur = zeros(1,xs) ;
    c_gra = zeros(1,xs) ;
    Point_Ecur = zeros(1,xs) ;
    Etot = 0 ;    
    for i = 1:xs 
        if(Tuv(i) == 1)
            acc = fix([xp(i) yp(i)]);
            Econ = ones(5,5) ;
            Ecur = ones(5,5) ;
            Egra = ones(5,5) ;
            min_m = ws+2 ;
            max_m = -(ws+2) ;
            min_n = ws+2 ;
            max_n = -(ws+2) ;
            n1 = ne1(i) ;
            n2 = ne2(i) ;

            Point_Ecur(i) = ((xp(n1) - 2*acc(1) + xp(n2))^2 + (yp(n1) - 2*acc(2) + yp(n2))^2) ;
            Point_Econ(i) = (d_mean-((xp(n1) - acc(1))^2 + (yp(n1) - acc(2))^2)^(1/2))^2 ;
            Point_Egra(i) =  Igrad(acc(1),acc(2)) ;

            for m = -ws:ws
                for n = -ws:ws
                    ac = acc + [m n] ;
                    ac = fix(ac) ;
                    if(ac(1) > 0 && ac(1) < rn && ac(2) > 0 && ac(2) < cn) 
                        Econ(m+ws+1,n+ws+1) = (d_mean-((xp(n1) - ac(1))^2 + (yp(n1) - ac(2))^2)^(1/2))^2 ;
                        Ecur(m+ws+1,n+ws+1) = ((xp(n1) - 2*ac(1) + xp(n2))^2 + (yp(n1) - 2*ac(2) + yp(n2))^2) ;
                        Egra(m+ws+1,n+ws+1) = Igrad(ac(1),ac(2)) ;
                        min_m = min([min_m m]) ;
                        max_m = max([max_m m]) ;
                        min_n = min([min_n n]) ;
                        max_n = max([max_n n]) ;
                    end
                end
            end
    %         min_m+ws+1:max_m+ws+1 ; min_n+ws+1:max_n+ws+1
            Econ = Econ(min_m+ws+1:max_m+ws+1,min_n+ws+1:max_n+ws+1) ;
            Ecur = Ecur(min_m+ws+1:max_m+ws+1,min_n+ws+1:max_n+ws+1) ;
            Egra = Egra(min_m+ws+1:max_m+ws+1,min_n+ws+1:max_n+ws+1) ;

            if(max(Econ(:)) ~= min(Econ(:)))
        %       Econ_Te{i} = (Econ)/max(Econ(:)) ;
                Econ_Te{i} = (Econ - min(Econ(:)))/(max(Econ(:)) - min(Econ(:))) ;
            else
                Econ_Te{i} = zeros(size(Econ)) ;
            end

            if(max(Ecur(:)) ~= min(Ecur(:)))
        %       Ecur_Te{i} = (Ecur)/max(Ecur(:)) ;
                Ecur_Te{i} = (Ecur - min(Ecur(:)))/(max(Ecur(:)) - min(Ecur(:))) ;
            else
                Ecur_Te{i} = zeros(size(Ecur)) ;
            end            
            if(max(Egra(:)) ~= min(Egra(:)))
                Egra_Te{i} = (Egra - min(Egra(:)))/(max(Egra(:)) - min(Egra(:))) ;
            else
                Egra_Te{i} = zeros(size(Egra)) ;
            end
        end
    end

    Point_Ecur = Point_Ecur / max(Point_Ecur(:)) ;
%    Point_Egra = Point_Egra / max(Point_Egra(:)) ;

    for i = 1:xs
        if(Tuv(i) == 1)
            n1 = ne1(i) ;
            n2 = ne2(i) ;
            c_con(i) = Con_We ; c_cur(i) = Cur_We ; c_gra(i) = Gra_We ;
                if((Point_Ecur(i) > Point_Ecur(n1)) & (Point_Ecur(i) > Point_Ecur(n2)) & (Point_Ecur(i) > cur_Tresh))
                    c_cur(i) = 0.2 ;
                end
        end
    end
    
    for i = 1:xs
        if(Tuv(i) == 1)
            Econ = Econ_Te{i} ;
            Ecur = Ecur_Te{i} ;
            Egra = Egra_Te{i} ;
            Eloc = c_con(i)*Econ + c_cur(i)*Ecur - c_gra(i)*Egra ;
            [nepx,nepy] = find(Eloc == min(Eloc(:))) ;
            xp(i) = nepx(1) + xp(i) - (ws+1) ;
            yp(i) = nepy(1) + yp(i) - (ws+1) ;
            Etot = Etot + min(Eloc(:)) ;
            if ( xp(i) < 1 | xp(i)> rn | yp(i) < 1 | yp(i) > cn )
                Tuv(i) = 0 ;
            end
        end
    end       
    New_Etot = Etot ;
end
m=1;
top_line=[2,length(lix)]
bottom_line=[2,length(lix)]
for m=1:length(lix)
     top_line(1,m)=lix(1,m);
     top_line(2,m)=liy(1,m);
     bottom_line(1,m)=lix(2,m);
     bottom_line(2,m)=liy(2,m);
     m+1;
 end
x_top=top_line(1,:);
y_top=top_line(2,:);

x_bottom=bottom_line(1,:);
y_bottom=bottom_line(2,:);

x_connect = [];
x_connect(length(x_bottom)+1) = x_bottom(1);
y_connect = [];
y_connect(length(y_bottom)+1) = y_bottom(1);

for index = 1:length(x_bottom)
    x_connect(index) = x_bottom(index);
    y_connect(index) = y_bottom(index);
end

figure, plot(x_connect, y_connect);

[img_y, img_x] = size(I);
snake=zeros(img_y,img_x);
%% 

% %draw in left breast pixel by pixel
% for a = 1:length(liy)                %a,b,n,m are just used as counters in the for loops - delete at end of section
%     e1 = liy;
%     for b = 1:length(e1.XData)              %get x and y data from cell array of ellipses and round so we can use them as indices
%         xe1(a,b) = e1.XData(b);
%         xe1(a,b) = round(xe1(a,b));
%         ye1(a,b) = e1.YData(b);
%         ye1(a,b) = round(ye1(a,b));
%         if xe1(a,b)<0.5
%             xe1(a,b)=1;
%         end
%         if ye1(a,b)<0.5
%             ye1(a,b)=1;
%         end
%     end
%     for d = 1:length(xe1)
%         snakes(ye1(a,d),xe1(a,d)) = 1;      %fill in 1's wherever there is a point in the ellipse
%     end
%     
% end
% [r,c]=size(I);
% for r=1:length(I)
%     a=find(r>lix && c>liy)
%     r+1;
% end
    
        
    

figure, imshow(I); line(lix,liy,'Color',[1 0 0],'linewidth',1);




%% Place Line on Top of Original Image from Kefset

ptID = input('Enter image name you want to open: ','s'); %Request patient image name
ptID = strcat(ptID,'.tif'); 
%I = imread(['C:\Users\smbailes\Documents\GitHub\loewlab\Segmentation MATLAB Code\Images\' ptID]); %open the image, keeping it in 16-bits
dir = uigetdir; 
I = imread([dir '/' ptID]); 

figure, imshow(I,[]) %to help decide if it should be cropped or not
title('Segmented Patient 10')
line(lix, liy,'Color',[1 0 0],'linewidth',1)

%% Crop Outside of Lines
%feed in lix and liy as points for which it can crop outside of
%look up MatLab cropping
cropimg = zeros(size(I));
for i = 1:length(x_connect)
    index_x = round(x_connect(i));
    index_y = round(y_connect(i));
    for k = 1:640
        for m = 1:480
            if k == index_x & m == index_y
                cropimg(m,k) = 1;
            end
        end
    end
end

figure, imshow(cropimg);
hold on;
fill(x_connect,y_connect,'w');
hold off;