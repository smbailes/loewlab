function CC = CONNCOMP(ImageMatrix, percent)

I = getMatrixOutliers(ImageMatrix);   % Remove Outliers
I_adj = I(find(I>0));       % Remove Zero Pixels
percent2 = percent / 100;
% % [b, edge] = histcounts(I_adj); % Get Image Histogram Data
% % Display Image Histogram w/ Percent Area
I_sort = sort(I_adj);       % Arrange Image Hist in Order Low -> High
percent_ind = round(percent2 * numel(I_sort));   % Find the index number for the User Input Percentage
percent_val = I_sort(end - percent_ind);        % Find Intensity for the Percentage Number
[overlay_r,overlay_c] = find(I >= percent_val); % Get Pixel locations above Percent Indicated

I_matr = I;
A = find(I_matr>=percent_val);
B = find(I_matr<percent_val);
I_matr(A) = true; 
I_matr(B) = false; 
CC = bwconncomp(I_matr);

end 