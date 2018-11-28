function centerline = symmetry_centerline(imageMatrix)
    [r,c] = size(imageMatrix); %Get Size of Current Image

    boundL = floor(c / 4); %Search Range for Midline is Middle 1/2 of the Image (Columnwise)
    boundU = ceil((3*c) / 4);

    colStart = imageMatrix(:,boundL); %Start Search Column from CurrentImage
    colEnd = imageMatrix(:,boundU); %End Search Column from Current Image

    holdVal = numel(colStart(find(colStart > 0))); %Holdval is the number of non-zero indices in the First column

    for i = boundL:boundU % Check only in the middle region of the image for the smallest non-zero column
        column = imageMatrix(:,i);

        check = numel(column(find(column > 0))); %Check is number of non-zero pixels in current column
        count(i) = check;

        if check < holdVal %hold is the column with the fewest nonzero pixels
            holdVal = check; %If current column has fewer nonero pixels, replace holdval with current column number
        end            

        xcent = floor(median(find(count == holdVal))); %If numerous columns have same number of minimal pixels, take middle column    
    end

    if ((xcent >= (boundU - 10)) || (xcent <=  (boundL + 10))) % if Center is found at borders of search range, reset image center to middle column
       xcent = round(c / 2); 
    end
    
    centerline = xcent;
end
