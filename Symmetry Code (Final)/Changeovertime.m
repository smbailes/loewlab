close all
clear, clc
answer = questdlg('ID Patient Type:','Patient Type','Patient','Volunteer','Patient');
if (strcmp(answer, 'Patient'))
    ptID = patientselect;    % Dialog Box for patient selection
    prompt = {'Enter User name:','Enter length of square in pixels:','Enter total number of pictures:'};
    dlgtitle = 'Input';
    defaultans = {'Jacob','20','14'};
    numlines = 1;
    answers = inputdlg(prompt,dlgtitle,numlines,defaultans);
    name = answers{1};
    location = (['C:\Users\' name '\Documents\GitHub\loewlab\Symmetry Code (Final)\Images\Patient Images\' ptID '\Cropped\']);
    
end
if (strcmp(answer, 'Volunteer'))
    vtID = volunteerselect;
    prompt = {'Enter User name:','Enter length of square in pixels:','Enter total number of pictures:'};
    dlgtitle = 'Input';
    defaultans = {'Jacob','20','14'};
    numlines = 1;
    answers = inputdlg(prompt,dlgtitle,numlines,defaultans);
    name = answers{1};
    location = (['C:\Users\' name '\Documents\GitHub\loewlab\Symmetry Code (Final)\Images\Volunteer Images\' vtID '\Cropped\']);
    
end
    
    %% Image Input
    numpics = str2double(answers{3}); % allocates number of pictures
    % Read images to cell matrix I_mat
    a=120; % set equal to the number of the first picture
    I_mat = cell(1,numpics);
    for i=1:numpics % set equal to total number of pictures being ran          
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
squareside = str2double(answers{2}); %converts ans to a number
for k = 1:numpics
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
for k = 1:numpics
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
stdv{i,j,k} = std2(imcrop(I_mat{k},square)); % takes the std2 of each block

    end
end
end
%% find good data via standard deviations
for i = 1:numrows
    for j = 1:numcols
        for k = 1:numpics
            if stdv{i,j,k} > 500 %possibly include an if loop to eliminate entire columns
                line{i,j,k} = '--'; %used in 
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
        for k = 1:numpics 
            if abs(stdv{i,j,k}) >500  % eliminates data that deviates too, mostly edge squares
                for d = 1:numpics % total number of pictures
                averages{i,j,d} = NaN;
                end   
            elseif averages{i,j,k} == 0
                for d = 1:numpics
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
        for k = 1:numpics
            if isempty(ylim_array{i,j,k}) % makes empty cells NaN
                ylim_array{i,j,k} = NaN;
            elseif ylim_array{i,j,k} <= 500 
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
%% graph averages. flat out destroys computers with very small squares
warning('off')
t = 2:numpics; % pictures start at t = 0, first pictue is lower than second. Strange
ypoints = cell(1,numpics); % preallocates y points. Might need to change
for i = 1:numrows
   figure(numpics + 1); 
    for j = 1:numcols
    ypoints = {};
        for k = 2:numpics % need to change
        ypoints{k} = averages{i,j,k};      
        end
      if isnan(ypoints{k}) == 0 %determines if the data is good for graphing       
        ypoints = cell2mat(ypoints); legend('show')
        coefficients = polyfit(t,ypoints,3); % creates the coefficients of the fitted curve. degree changes
        newypoints = polyval(coefficients,t); % creates new y points that are smoooth
        if sqrt(numrows) - floor(sqrt(numrows))  <0.5
            k = 1;
        else 
            k = 0;
        end
        subplot(ceil(sqrt(numrows)),ceil(sqrt(numrows)),i) % creates subplot
        r = rand;, g = rand;, b = rand; %sets random rgb values to make lots of colots
        plot(t,newypoints,'-','Color',[r g b],'DisplayName',num2str(j)), hold on %can add real data points with t,y,'+'
        g =  plot(t,ypoints,'o','Color',[r g b]); %Adds true data points to graph
        gAnno = get(g, 'Annotation'); %following three lines make the true data points not appear in legend
        gLegend = get(gAnno, 'LegendInformation');
        set(gLegend,'IconDisplayStyle','off');
        title(['row ' num2str(i)])
        xlabel('time')
        ylabel('pixel value')
        ylim([ymin,ymax]); %specify y limits
        xlim([2,numpics]); % need to change
      end
    end
end
%% Find the average slope of each square (ohgod)
for i = 1:numrows
    for j = 1:numcols
        for k = 1:numpics -1
            change{i,j,k}= averages{i,j,k+1} - averages{i,j,k};
        end
    end
end
avechange = cell(numrows,numcols)
for i = 1:numrows
    for j =1:numcols
        avechange{i,j} = nanmean(cell2mat(change(i,j,:)))
    end
end
%% Find the total change of each square
totsquarechange = cell(numrows,numcols);
for i = 1:numrows
    for j = 1:numcols
        totsquarechange{i,j} = averages{i,j,numpics} - averages{i,j,2};
    end
end
   
%% graph Standard deviations

% t = [0:numpics - 1];
% for i = 1:numrows
%     figure
%     for j = 1:numcols
%     xpoints = {};
%         for k = 1:numpics
%    
%         xpoints{k} = stdv{i,j,k};
%         end
%     xpoints = cell2mat(xpoints);
%     plot(t,xpoints), hold on
%     legend({'col 1','col 2', 'col 3', 'col 4', 'col 5', 'col 6', 'col 7', 'col 8', 'col 9', 'col 10'})
%     title(['Standard Deviation: row ' num2str(i)])
%     xlabel('time')
%     ylabel('pixel value')
%     end
% 
%  end
%% finding change from start to finish. Isn't helpful with large squares
for i = 1:numrows
    for j = 1:numcols
        totchange{i,j} = (averages{i,j,numpics} - averages{i,j,2}); %creats a cell array of the change in T for each square
        if totchange{i,j} > 0
            totchange{i,j} = NaN;
        end
    end
end

[xpoints,ypoints] = meshgrid(1:numcols,1:numrows);
zpoints = cell2mat(totchange);
figure(numpics + 2), surface(xpoints,ypoints,zpoints)
view(-37,64);
xlabel('Column')
ylabel('Row')
zlabel('Change in pixel value')
title('Relative Change in Temperature Of Sections Of The Breasts')
axis ij % makes axis match the figures
 %% find the averages of each breast
prompt = ('Enter middle column');
dlgtitle = ('Middle column');
numlines = 1;
figure(numpics)
midcol = inputdlg(prompt,dlgtitle,numlines);

for i = 1:numrows % Takes data from right breast
    for j = 1:str2double(midcol{1})
        for k = 2:numpics % will have start at 1 
            Rbreastmean{i,j,k} = averages{i,j,k};
            if Rbreastmean{i,j,k} == 0
                Rbreastmean{i,j,k} = NaN;
            end
        end
    end 
end
totRbreastmean = cell(1,numpics -1); % allocate arrays
Rstdv = cell(1,numpics-1);
for k =2:numpics % will have to get rid of -1. Finds averages of Rbreast and stdv 
    rightmean = cell2mat(Rbreastmean(:,:,k));
    rightmean = nanmean(rightmean);
    totRbreastmean{k} = nanmean(rightmean);
    Rightstdv = cell2mat(Rbreastmean(:,:,k));
    Rightstdv = nanstd(Rightstdv);
    Rstdv{k} = nanstd(Rightstdv);
end
for i = 1:numrows %gets the data for the left breast
    for j = str2double(midcol{1}):numcols
        for k = 2:numpics % will have start at 1 
            Lbreastmean{i,j,k} = averages{i,j,k};
            if Lbreastmean{i,j,k} == 0
                Lbreastmean{i,j,k} = NaN;
            end
        end
    end 
end
for i = 1:numrows %Makes all zero cells = NaN
    for j = 1:str2double(midcol{1})-1
        for k = 2:numpics % will have start at 1 
            Lbreastmean{i,j,k} = NaN;
            
        end
    end 
end
totLbreastmean = cell(1,numpics -1);
Lstdv = cell(1,numpics-1);
for k =2:numpics % will have to get rid of -1. calculates mean and stdv
    leftmean = cell2mat(Lbreastmean(:,:,k));
    leftmean = nanmean(leftmean);
    totLbreastmean{k} = nanmean(leftmean);
    Leftstdv = cell2mat(Lbreastmean(:,:,k));
    Leftstdv = nanstd(Leftstdv);
    Lstdv{k} = nanstd(Leftstdv);
end

figure(numpics + 3) % get everything into a matrix
t = 2:numpics;
totRbreastmean = cell2mat(totRbreastmean);
totLbreastmean = cell2mat(totLbreastmean);
Lstdv = cell2mat(Lstdv);
Rstdv = cell2mat(Rstdv);


Rplot = errorbar(t,totRbreastmean,Rstdv,'r');, hold on % plot everything
Lplot = errorbar(t,totLbreastmean,Lstdv,'b');

title('Left vs Right Change Over Time')
legend('Right','Left')
xlabel('Time (min)')
ylabel('Pixel Value')
a = Rplot.Color;
Rplot.Color = 'r'
b = Lplot.Color;
Lplot.Color = 'b';

%% find average slope for each breast
for k = 1:numpics-2 % will have to make 1
    Lbreastchange{k} =  totLbreastmean(k+1) - totLbreastmean(k);
end
aveLbreastchange = mean(cell2mat(Lbreastchange));
for k = 1:numpics-2 % will have to make 1
    Rbreastchange{k} =  totRbreastmean(k+1) - totRbreastmean(k);
end
aveRbreastchange = mean(cell2mat(Rbreastchange));
%% Find total change for each breast
totLbreastchange = totLbreastmean(numpics-1) - totLbreastmean(1);
totLbreastchange = totRbreastmean(numpics-1) - totRbreastmean(1);

