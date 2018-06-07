path = uigetdir;

path = strcat(path, '\0120.tif');
I = imread(path);
Igpu = I;

[r c] = find(Igpu);


for a = 1:length(r)
    region(a,1) = c(a); %X-coordinate
    region(a,2) = r(a); %Y-coordinate
    region(a,3) = I(r(a), c(a)); %value at that x,y coordinate
end  



% for r = 1:rows
%     for c = 1:cols
%         I_xyz(r, c, i) = Igpu(r, c);
%         i = i+1;
%     end
% end


