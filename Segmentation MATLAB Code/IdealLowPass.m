
% fc is the circular cutoff frequency which is normalized to [0 1], that is, 
% the highest radian frequency \pi of digital signals is mapped to 1.

ptID = input('Enter image name you want to open: ','s'); % ENTER PATIENT IMAGE NUMBER
%ex. 0000, 0120, 0240, 0360, 0480...1800 or 1799
ptID = strcat(ptID,'.tif'); % ptID labelled as TIF file (ex. 0000.tif)
dir = uigetdir; % select patient number's folder from folder with 8 bit images
im0 = imread([dir '/' ptID]); % retrieves image from selected folder 

fc = 0.75;

figure(1), imshow(im0, [])

[ir,ic,iz] = size(im0); 
hr = (ir-1)/2; 
hc = (ic-1)/2; 
[x, y] = meshgrid(-hc:hc, -hr:hr);

mg = sqrt((x/hc).^2 + (y/hr).^2); 
lp = double(mg <= fc);

IM = fftshift(fft2(double(im0))); 
IP = zeros(size(IM)); 
for z = 1:iz 
IP(:,:,z) = IM(:,:,z) .* lp; 
end 
im = abs(ifft2(ifftshift(IP)));

figure(2), imshow(im, [])

