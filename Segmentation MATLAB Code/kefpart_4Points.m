% Point system code.
% The idea is to look pixel by pixel to determine which 

final_im=zeros(img_y,img_x);
counter=1;
while counter<5
bound_canny=bwmorph(boundededges,'bridge');
canny_bounded=bwmorph(bound_canny,'thicken');
counter=counter+1;
end

for aa = 1:img_y
    for bb = 1:img_x
        if newI(aa,bb)~=0
            final_im(aa,bb)=final_im(aa,bb)+1;
        end
        if lowerellipse(aa,bb)==1
            final_im(aa,bb)=final_im(aa,bb)+2;
        end
        if canny_bounded(aa,bb)==1
            final_im(aa,bb)=final_im(aa,bb)+1;
        end
        if circs(aa,bb)==1
            final_im(aa,bb)=final_im(aa,bb)+1;
        end
        if I_sobel(aa,bb)==1
            final_im(aa,bb)=final_im(aa,bb)+2;
        end
    end
end

figure,imshow(final_im,[])
title('Combined Point Systems')


%Only the overlapped

finfin = zeros(img_y,img_x);
for cc=1:img_y
    for dd=1:img_x
        if final_im(cc,dd)>2
            finfin(cc,dd)=1;
        end
    end
end

clear aa bb cc dd
% 
figure, imshow(I,[]), title('Breast Pre-Cleaning')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', finfin) 

know=input('More cleaning required? [y/n]: ','s')';
if know == 'y'
    kefpart_6Clean
else
    zm_SB7
end

