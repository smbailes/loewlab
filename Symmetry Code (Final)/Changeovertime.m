close all
clear, clc
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
     
%% create grid over image and find averages
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
            if stdv{i,j,k} > 500 %possibly include an if loop to eliminate entire columns
                line{i,j,k} = '--';
               % averages{i,j,k} = NaN; %eliminates bad data
            else 
                line{i,j,k} = '-';
                gooddata{i,j,k} = averages{i,j,k}; %separates good data into separate array
                
            end
        end
    end
end
%% find good data for 

for i = 1:numrows
    for j = 1:numcols
        for k = 1:14 %number of total pictures - 1
            if abs(averages{i,j,k} - averages{i,j,k+1}) >500  %eliminates data that deviates too much at any one point
                for d = 1:15 % total number of pictures
                averages{i,j,d} = NaN;
                end   
            elseif averages{i,j,k} == 0
                for d = 1:15
                    averages{i,j,k} = NaN;
                end
            end
        end
    end
end

%% determining y-value limits
ylim_array = gooddata;
for i = 1:numrows
    for j = 1:numcols
        for k = 1:15
            if isempty(ylim_array{i,j,k});
                ylim_array{i,j,k} = NaN;
            elseif ylim_array{i,j,k} <= 500;
                ylim_array{i,j,k} = NaN;
            elseif ylim_array{i,j,k} == 0
                for d = 1:15
                    ylim_array{i,j,d} = NaN;
                end
            end
        end
    end
end
ymin = min(cell2mat(ylim_array));
ymin = min(ymin);
ymin = min(ymin);
ymax = max(cell2mat(ylim_array));
ymax = max(ymax);
ymax = max(ymax);
%% graph averages
warning('off')
t = 0:14;
counter = 1;
jcounter = 0;
%colcount = zeros(1,numcols)
for i = 1:numrows
    figure; 
    for j = 1:numcols
    ypoints = {};
        for k = 1:15
        ypoints{k} = averages{i,j,k};      
        end
      if isnan(ypoints{k}) == 0 %determines if the data is good for graphing
        % xpoints = {NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN};
        
        ypoints = cell2mat(ypoints); legend('show')
        coefficients = polyfit(t,ypoints,3); % creates the coefficients of the fitted curve
        newypoints = polyval(coefficients,t); % creates new y points that are smoooth
        plot(t,newypoints,'-','DisplayName',num2str(j)), hold on %line{i,j,k} can be used in place of *
        title(['change of regions over time: row ' num2str(i)])
        xlabel('time')
        ylabel('pixel value')
        ylim([ymin,ymax]); %specify y limits
      end
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