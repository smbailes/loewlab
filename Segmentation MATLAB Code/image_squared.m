close all;

for a = 1:img_y
    for b = 1:img_x
        I2(a,b) = (I(a,b))*0.01;
    end
end

figure(1)
imshow(I,[])

figure(2)
imshow(I2,[])

edgecanny1 = edge(I,'canny');
edgecanny1=bwareaopen(edgecanny1,10); %removes very small edge lines

figure,imshow(edgecanny1)
title('Canny edges 1');

edgecanny2 = edge(I2,'canny');
edgecanny2=bwareaopen(edgecanny2,10); %removes very small edge lines

figure,imshow(edgecanny2)
title('Canny edges 2');