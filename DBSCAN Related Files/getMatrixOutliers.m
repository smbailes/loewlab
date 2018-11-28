function newmat = getMatrixOutliers(matrix) 
%% Using Quartile Method - ONLY REMOVING POS. END OUTLIERS (> 1.5*IQR)  
%% REPLACING OUTLIERS WITH MEDIAN OF SURROUNDING PIXELS
    [nr, nc] = size(matrix);
    A = reshape(matrix,[1, numel(matrix)]); %make matrix into single vector
    B = sort(A); %Arrange vector of image intensities from smallest to largest
    
    ind75 = round(0.75*numel(B)); %Index of 75th percentile pixel
    ind25 = round(0.25*numel(B)); %Index of 25th percentile pixel
        
    quartile75 = B(ind75); %Pixel Intensity at 75th percentile
    quartile25 = B(ind25); %Pixel intensity at 25th percentile
    iqr = quartile75 - quartile25; %Interquartile Range of Intensity
    
    
    % OUTLIER IS DEFINED AS ANY PIXEL WITH INTENSITY ABOVE 1.5x THE IQR
    % INCREASE VALUE FOR HIGHER THRESHOLD (LESS PIXELS REMOVED)
    % DECREASE VALUE FOR LOWER THRESHOLD (MORE PIXELS REMOVED)
    D = 1.5; %THRESHOLD VALUE
    [y_pos,x_pos] = find(matrix > (quartile75 + (D*iqr))); %CHANGE HERE FOR DIFFERENT OUTLIER STANDARDS.  
    newmat = matrix;
   
    for i = 1:(length(x_pos))
        if (x_pos(i) == 1)&&(y_pos(i) == 1) %TOP LEFT CORNER
            newmat(y_pos(i),x_pos(i)) = median([matrix(y_pos(i),x_pos(i)+1)  %Replace with Median of surrouding pixels
                                        matrix(y_pos(i),x_pos(i)+1) 
                                        matrix(y_pos(i)+1,x_pos(i)+1)]);
        
        elseif (x_pos(i) == nc) && (y_pos(i) == 1) %TOP RIGHT CORNER
            newmat(y_pos(i),x_pos(i)) = median([matrix(y_pos(i),x_pos(i)-1)  %Replace with Median of surrouding pixels
                                        matrix(y_pos(i)+1,x_pos(i)-1) 
                                        matrix(y_pos(i)+1,x_pos(i))]);
        
        elseif (x_pos(i) == 1) && (y_pos(i) == nr) %BOTTOM LEFT CORNER
            newmat(y_pos(i),x_pos(i)) = median([matrix(y_pos(i)-1,x_pos(i))  %Replace with Median of surrouding pixels
                                        matrix(y_pos(i)-1,x_pos(i)+1) 
                                        matrix(y_pos(i),x_pos(i)+1)]);
        
        elseif (x_pos(i) == nc) && (y_pos(i) == nr) %BOTTOM RIGHT CORNER
            newmat(y_pos(i),x_pos(i)) = median([matrix(y_pos(i)-1,x_pos(i)-1)  %Replace with Median of surrouding pixels
                                        matrix(y_pos(i)-1,x_pos(i)) 
                                        matrix(y_pos(i),x_pos(i)-1)]);
        
        elseif (x_pos(i) == 1) && (y_pos(i) > 1 && y_pos(i) < nr) %ALONG FIRST COLUMN
            newmat(y_pos(i),x_pos(i)) = median([matrix(y_pos(i)-1,x_pos(i))  %Replace with Median of surrouding pixels
                                        matrix(y_pos(i)-1,x_pos(i)+1)
                                        matrix(y_pos(i),x_pos(i)+1) 
                                        matrix(y_pos(i)+1,x_pos(i)) 
                                        matrix(y_pos(i)+1,x_pos(i)+1)]);
        
        elseif (x_pos(i) == nc) && (y_pos(i) > 1 && y_pos(i) < nr) %ALONG LAST COLUMN
            newmat(y_pos(i),x_pos(i)) = median([matrix(y_pos(i)-1,x_pos(i)-1)  %Replace with Median of surrouding pixels
                                        matrix(y_pos(i)-1,x_pos(i)) 
                                        matrix(y_pos(i),x_pos(i)-1) 
                                        matrix(y_pos(i)+1,x_pos(i)-1)  
                                        matrix(y_pos(i)+1,x_pos(i))]);
        
        elseif (y_pos(i) == 1) && (x_pos(i) > 1 && x_pos(i)< nc) %ALONG FIRST ROW
            newmat(y_pos(i),x_pos(i)) = median([matrix(y_pos(i),x_pos(i)-1)  %Replace with Median of surrouding pixels
                                        matrix(y_pos(i),x_pos(i)+1) 
                                        matrix(y_pos(i)+1,x_pos(i)-1) 
                                        matrix(y_pos(i)+1,x_pos(i)) 
                                        matrix(y_pos(i)+1,x_pos(i)+1)]);
        
        elseif (y_pos(i) == nr) && (x_pos(i) > 1 && x_pos(i)< nc) %ALONG LAST ROW
            newmat(y_pos(i),x_pos(i)) = median([matrix(y_pos(i)-1,x_pos(i)-1)  %Replace with Median of surrouding pixels
                                        matrix(y_pos(i)-1,x_pos(i)) 
                                        matrix(y_pos(i)-1,x_pos(i)+1) 
                                        matrix(y_pos(i),x_pos(i)-1) 
                                        matrix(y_pos(i),x_pos(i)+1)]);
        
        elseif (x_pos(i) > 1 && x_pos(i)< nc) && (y_pos(i) > 1 && y_pos(i) < nr) %ANY LOCATION MATRIX INTERIOR
            newmat(y_pos(i),x_pos(i)) = median([matrix(y_pos(i)-1,x_pos(i)-1)  %Replace with Median of surrouding pixels
                                        matrix(y_pos(i)-1,x_pos(i)) 
                                        matrix(y_pos(i)-1,x_pos(i)+1) 
                                        matrix(y_pos(i),x_pos(i)-1) 
                                        matrix(y_pos(i),x_pos(i)+1) 
                                        matrix(y_pos(i)+1,x_pos(i)-1) 
                                        matrix(y_pos(i)+1,x_pos(i)) 
                                        matrix(y_pos(i)+1,x_pos(i)+1)]); %Replace Outlier with mean of surrounding pixels
        end
    
    end
    
%     A(pos_outliers(end)) = A(pos_outliers(end) - 1);
    
%     for i = 1:length(neg_outliers)
%             A(neg_outliers(i)) = mean([A(neg_outliers(i)-1), A(pos_outliers(i)+1)]); %Replace Outlier with mean of surrounding pixels
%     end
    
end