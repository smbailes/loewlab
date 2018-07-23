clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 15;

% Read in a standard MATLAB gray scale demo image.
%folder = fullfile(matlabroot, '\toolbox\images\imdemos');
baseFileName = '0000 - P12.tif';
% baseFileName = 'eight.tif';
% Get the full filename, with path prepended.
%fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
% if ~exist(fullFileName, 'file')
% 	% File doesn't exist -- didn't find it there.  Check the search path for it.
% 	fullFileName = baseFileName; % No path this time.
% 	if ~exist(fullFileName, 'file')
% 		% Still didn't find it.  Alert user.
% 		errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
% 		uiwait(warndlg(errorMessage));
% 		return;
% 	end
% end
grayImage = imread(baseFileName);
% Get the dimensions of the image.  
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(grayImage);
if numberOfColorBands > 1
	% It's not really gray scale like we expected - it's color.
	% Convert it to gray scale by taking only the green channel.
	grayImage = grayImage(:, :, 2); % Take green channel.
end
% Display the original gray scale image.
subplot(2, 2, 1);
imshow(grayImage, []);
axis on;
title('Original Grayscale Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 
drawnow;

% Let's compute and display the histogram.
[pixelCount, grayLevels] = imhist(grayImage);
subplot(2, 2, 2); 
bar(grayLevels, pixelCount);
grid on;
title('Histogram of original image', 'FontSize', fontSize);
xlim([0 grayLevels(end)]); % Scale x axis manually.

% Construct an initial mask that is the convex hull of everything.
subplot(2, 2, 2);
% mask = grayImage < 175; % Use for eight.tif
mask = grayImage > 80; % Use for coins.png
mask = bwconvhull(mask, 'Union');
imshow(mask);
axis on;
title('Initial mask', 'FontSize', fontSize);
drawnow;

% Now find the improved outer boundary using active contours.
bw = activecontour(grayImage, mask, 10000);
subplot(2, 2, 3);
imshow(bw);
axis on;
title('Outer Boundary Mask', 'FontSize', fontSize);
drawnow;

% Display the original gray scale image in the lower left
% So we can display boundaries over it.
subplot(2, 2, 4);
imshow(grayImage, []);
axis on;
hold on;

% Display the initial contour on the original image in blue.
% Display the final contour on the original image in red.
contour(mask,[0 0],'b', 'LineWidth', 4); 
contour(bw,[0 0],'r', 'LineWidth', 4); 
title('Image with initial and final contours overlaid', 'FontSize', fontSize);

msgbox('Done with activecontour demo');