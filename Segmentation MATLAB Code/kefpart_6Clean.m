finfin(1:round(img_y/4),:)=0;


conti=1;
while conti<9
    Connect3d=bwconncomp(finfin);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [smallest,idx] = min(numPixels);
    finfin(CC.PixelIdxList{idx}) = 0;
    conti=conti+1;
end




%This part of the code is to clean up what has been created using the point
%system in the previous parts. 
remov = bwareaopen(finfin, 5);

closer = bwmorph(remov,'bridge');
%figure,imshow(closer)

countit = 1;
while countit < 7
    closer = bwmorph(closer,'bridge');
    countit=countit+1;
end



finfin = bwmorph(closer,'bridge');
finfin=bwmorph(finfin,'bridge');
finfin = bwareaopen(finfin,50);
figure, imshow(I,[]), title('Breast Pre-Cleaning')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
set(displ, 'AlphaData', finfin) 


% edgy=edge(closer,'canny');
% %figure,imshow(edgy)

% 
% if in == 's'
%     remover = bwareaopen(edgy, 5);
% else
%     remover = bwareaopen(edgy, 5);
% end

% [e,f]=size(I);
% tryit=zeros(e,f);
% for fir=1:f
%     y=find(remover(:,fir)==1);
%     if ~isempty(y)
%     tryit(y(1),fir)=1;
%     end
% end

finfin=bwmorph(finfin,'bridge');


figure, imshow(I,[]), title('Final')
% blue on top on figure
blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
hold on 
displ = imshow(blue); 
hold off 
% Use our diff1 as the AlphaData for the solid red image. 
set(displ, 'AlphaData', finfin) 

zm_7_logedges

% kefpart_7finishconnections

% 
% sra=connecttt;
% %to get rid of extraneous higher pixels:
% [xxx,yyy]=find(connecttt(:,tmp1:tmp2)==1);
% highest=min(xxx);
% sra(1:highest,:)=0;
% 
% 
% andnow = bwareaopen(sra, 3);
% 
% 
% figure, imshow(I,[]), title('Final')
% % blue on top on figure
% blue = cat(3, zeros(size(I)), zeros(size(I)), ones(size(I))); %blue has RGB value 0 0 1
% hold on 
% displ = imshow(blue); 
% hold off 
% % Use our diff1 as the AlphaData for the solid red image. 
% set(displ, 'AlphaData', andnow) 