%% Component Isolation
% Finds mid line of image and isolates connected components to be used to
% crop image without taking random small lines on sides into account

% Last Updated by Zainab Mahmood on 2/5/18




%% Part 1: Split image in half and find connected components

% find connected components
CI_log = bwconncomp(bibi);

compi_x = cell(1,CI_log.NumObjects);     % create cell arrays for x and y coordinates for each conn comp
compi_y = cell(1,CI_log.NumObjects);
for n = 1:CI_log.NumObjects
    [compi_y{n}, compi_x{n}] = ind2sub(size(I),CI_log.PixelIdxList{n});               % turn pixel indices into x,y coordinates to graph later
end

% graph components separately
conncompimat = cell(1,CI_log.NumObjects);
for n = 1:CI_log.NumObjects
    conncompimat{n} = zeros(img_y,img_x);
    for a = 1:length(compi_x{n})                                % put each component in a matrix, which is part of a cell array
        conncompimat{n}(compi_y{n}(a),compi_x{n}(a)) = 1;         % each cell is a different component
    end
    
    % graph components
    figure;
    imshow(I,[]);
    mytitle = sprintf('Connected Component # %i',n);
    title(mytitle);
    red = cat(3, ones(size(I)), zeros(size(I)), zeros(size(I))); %red is 100
    hold on
    h7 = imshow(red);
    set(h7, 'AlphaData', conncompimat{n});
    hold off
end

% split image in half
figure;
imshow(I,[]);
hold on;
mid_col = zeros(img_y,img_x);
mid_col(:,img_x/2) = 1;
h6 = imshow(cyan);
set(h6, 'AlphaData', mid_col);
hold off

%% Part 2: Find and plot top pixel of each column




