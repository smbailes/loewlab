keepingtheold=logfin; %this way I can keep the original to work with so I don't have to go through the whole process.

findthin = find(bottoms3>0);
[thiny, thinx] = ind2sub(size(I),findthin);
[P,indP]=min(findthin); %finds the first indece to start with
[xout, yout] = points2contour(thinx,thiny,2,'ccw');
figure;
plot(xout,yout); %Plotting the xout and yout, not over the image.
title('Closed contour');
set(gca,'YDir','reverse');

% plot points xout,yout on image
contourmap = zeros(img_y,img_x);
for alpha = 1:length(xout)
    contourmap(yout(alpha),xout(alpha))= 2^16;
end

figure;
imshow(I,[]);
hold on;
comet(xout,yout);
axis([0 640 0 480]);
title('Xout, Yout');


% vertices = horzcat(yout',xout');
% curvature = LineCurvature2D(vertices);
% figure;
% plot(curvature);
% title('curvature');
% 
% 
% for fig=1:length(curvature)
%     if curvature(fig)==0 %takes out the pixels where the line is equal to zero.
%         logfin(vertices(fig,1),vertices(fig,2))=0;
%     end
% end



figure;
imshow(I,[]);
title('getting rid of 0 curvatures');
green = cat(3, ones(size(I)), zeros(size(I)), ones(size(I)));
hold on
h4 = imshow(green);
set(h4, 'AlphaData', logfin);
hold off



afterp2c=zeros(img_y,img_x); %attempting to put xout and yout into a binary image
for ty=1:length(xout)
    afterp2c(yout(ty),xout(ty))=2^16;
end


figure;
imshow(I,[])
hold on
plot(xout,yout,'b-')
set(gca,'YDir','reverse');
title('Plotted over the Image') %plotting xout and yout over the image

logaft=logical(afterp2c);
biggest = bwareafilt(logaft,1,'largest');


figure;
imshow(I,[]);
title('Implementing Xout and Yout');
green = cat(3, ones(size(I)), zeros(size(I)), ones(size(I)));
hold on
h4 = imshow(green);
set(h4, 'AlphaData', afterp2c);
hold off