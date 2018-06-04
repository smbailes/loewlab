close all
answer = questdlg('ID Patient Type:','Patient Type','Patient','Volunteer','Patient');
if (strcmp(answer, 'Patient'))
    ptID = patientselect;    % Dialog Box for patient selection
    %% User Selection (CHANGE THIS TO NEW USERS / DIRECTORY REF SYSTEM)
    user = userselect;          % Dialog Box for user selection
    a = strcmp(user,'Sydney');   % Compare user input to 'Sydney"
    b = strcmp(user,'Samhita');
    c = strcmp(user,'Jacob');
    if (a == 1 && b == 0)                 % If Sydney was selected
        location = (['C:\Users\smbailes\Documents\GitHub\loewlab\Symmetry Code (Final)\Images\' ptID '\Cropped Image\']);
    elseif (a == 0 && b == 1)  % If Samhita was selected
         location = (['C:\Users\samhitamurthy\GitHub\loewlab\Symmetry Code (Final)\Images\' ptID '\Cropped Image\']);
    elseif (a == 0 && b == 0 && c == 1) % If Jacob was selected
        location = (['C:\Users\Jacob\Documents\GitHub\loewlab\Symmetry Code (Final)\Images\' ptID '\Cropped Image\']);
    end
end
    %% Image Input
    
    % Read 15 images to cell matrix I_mat
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end
    
    I1 = I_mat{1};              % Display first image
    I = getMatrixOutliers(I1);  % Remove outliers
    I_adj1 = I1(find(I1>0));    % Remove zero pixels
     
%% create grid over image
prompt = ('Enter size of one box on grid in pixels'); % make dialog box
dlgtitle = ('input');
num_lines = 1;
defaultans = {'50'};
ans = inputdlg(prompt,dlgtitle,num_lines,defaultans);
squareside = str2num(cell2mat(ans)); %converts ans to a number
for k = 1:15
I_mat{k}(squareside:squareside:end,:,:) = 0;% converts every nth row to black
I_mat{k}(:,squareside:squareside:end,:) = 0;% converts every nth column to black

[r,c] = size(I_mat{k});
numrows = floor(r/squareside); %calculates the number of full rows
numcols = floor(c/squareside); %calculates the number of full columns


for j = 1:1:numcols
    row = squareside*(j-1)+1;
    for i = 1:1:numrows
        col = squareside*(i-1)+1 ;
square = [row, col, squareside-2, squareside-2]; %  creates the square to be averaged
averages{i,j,k} = mean2(imcrop(I_mat{k},square)); % takes the average of each block

    end
end
figure, imshow(I_mat{k},[min(I_adj1) max(I_adj1)]) % displays each image at each minute
end
%% standard Deviation
for k = 1:15
I_mat{k}(squareside:squareside:end,:,:) = 0;% converts every nth row to black
I_mat{k}(:,squareside:squareside:end,:) = 0;% converts every nth column to black

[r,c] = size(I_mat{k});
numrows = floor(r/squareside); %calculates the number of full rows
numcols = floor(c/squareside); %calculates the number of full columns


for j = 1:1:numcols
    row = squareside*(j-1)+1;
    for i = 1:1:numrows
        col = squareside*(i-1)+1 ;
square = [row, col, squareside-2, squareside-2]; %  creates the square to be averaged
stdv{i,j,k} = std2(imcrop(I_mat{k},square)); % takes the average of each block

    end
end
end
%% find good data via standard deviations
for i = 1:numrows
    for j = 1:numcols
        for k = 1:15
            if stdv{i,j,k} > 250
                line{i,j,k} = '--';
            else 
                line{i,j,k} = '-';
            end
        end
    end
end
%% graph averages
warning('off')
t = [0:14];
for i = 1:numrows
    figure
    for j = 1:numcols
    xpoints = {};
        for k = 1:15
   
        xpoints{k} = averages{i,j,k};
        end
    xpoints = cell2mat(xpoints);
    plot(t,xpoints,line{i,j,k}), hold on
    legend({'col 1','col 2', 'col 3', 'col 4', 'col 5', 'col 6', 'col 7', 'col 8', 'col 9', 'col 10'})
    title(['change of regions over time: row ' num2str(i)])
    xlabel('time')
    ylabel('pixel value')
    end

end
%% graph Standard deviations
%{
t = [0:14];
for i = 1:numrows
    figure
    for j = 1:numcols
    xpoints = {};
        for k = 1:15
   
        xpoints{k} = stdv{i,j,k};
        end
    xpoints = cell2mat(xpoints);
    plot(t,xpoints), hold on
    legend({'col 1','col 2', 'col 3', 'col 4', 'col 5', 'col 6', 'col 7', 'col 8', 'col 9', 'col 10'})
    title(['Standard Deviation: row ' num2str(i)])
    xlabel('time')
    ylabel('pixel value')
    end

end
%}