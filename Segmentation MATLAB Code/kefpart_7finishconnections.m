CC2 = bwconncomp(connected);
newboundaries4 = connected;

%Idea is to connect the last pixel from conncomp 1 with the first pixel on
%conncomp 2. The first and last pixels are determined by the coloums. Last
%pixel is the pixel with the max column index, and first pixel is with the
%min column index on the connected component. Remember that the columns
%indeces are numbered from left to right of the image.

for n = 1:CC2.NumObjects - 1 %the number of lines in the middle region of the patient
    %Store all row and col values of component n and the component after in
    % x1,y1, x2, y2
    [x1, y1] = ind2sub(size(newboundaries4),CC2.PixelIdxList{n});
    [x2, y2] = ind2sub(size(newboundaries4),CC2.PixelIdxList{n+1});
    
    [yy1, ind] = max(y1); %find the max col in component n
    xx1 = x1(ind); % The corrosponding row value for max col
    
    [yy2, ind] = min(y2); %find the min col in component n+1
    xx2 = x2(ind); % The corrosponding row value for min col
    
    %Draw a line between the two points (xx1,yy1) and (xx2,yy2) and insert
    %it in newboundaries4
    shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'BorderColor', 'White');
    newboundaries4 = step(shapeInserter, newboundaries4, uint16([yy1 xx1 yy2 xx2]));
    %figure, imshow(newboundaries4), title('After step shapeinserter');
    
end
clear xx2 xx1 yy1 yy2 y1 y2 x1 x2;


%newboundaries4=bwmorph(newboundaries4,'skel',Inf);


figure, imshow(I,[]), title('Connecting')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
% Use our diff1 as the AlphaData for the solid red image. 
set(displ, 'AlphaData', newboundaries4)


