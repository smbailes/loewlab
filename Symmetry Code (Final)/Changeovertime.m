close all
clc, clearvars -except normbreastmean normstdv tottumbreastmean tumstdv v
v = v+1
answer = questdlg('ID Patient Type:','Patient Type','Patient','Volunteer','Patient');
if (strcmp(answer, 'Patient'))
    ptID = patientselect;    % Dialog Box for patient selection
    prompt = {'Enter User name:','Enter length of square in pixels:','Enter total number of pictures:'};
    dlgtitle = 'Input';
    defaultans = {'Jacob','16','15'};
    numlines = 1;
    answers = inputdlg(prompt,dlgtitle,numlines,defaultans);
    name = answers{1};
    location = (['C:\Users\' name '\Documents\GitHub\loewlab\Symmetry Code (Final)\Images\Patient Images\' ptID '\Cropped\']);
    
end
if (strcmp(answer, 'Volunteer'))
    vtID = volunteerselect;
    prompt = {'Enter User name:','Enter length of square in pixels:','Enter total number of pictures:'};
    dlgtitle = 'Input';
    defaultans = {'Jacob','16','14'};
    numlines = 1;
    answers = inputdlg(prompt,dlgtitle,numlines,defaultans);
    name = answers{1};
    location = (['C:\Users\' name '\Documents\GitHub\loewlab\Symmetry Code (Final)\Images\Volunteer Images\' vtID '\Cropped\']);
    
end
    
    %% Image Input
    numpics = str2double(answers{3}); % allocates number of pictures
    % Read images to cell matrix I_mat
    a=0000; % set equal to the number of the first picture
    I_mat = cell(1,numpics);
    for i=1:numpics % set equal to total number of pictures being ran          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end
    
    I1 = I_mat{1};              % Display first image
    I1 = getMatrixOutliers(I1);  % Remove outliers
    I_adj1 = I1(find(I1>0));   % Remove zero pixels
    I = cell(1,numpics);
    for k = 1:numpics
        I{k} = I_mat{k};
    end
     
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
axis on
xticks([squareside/2:squareside:c]) %adds axes to images
xticklabels([1:1:numcols])
yticks([squareside/2:squareside:c])
yticklabels([1:1:numrows])
set(gca,'XaxisLocation','top')
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
t = 0:numpics-1; % pictures start at t = 0, first pictue is lower than second. Strange
ypoints = cell(1,numpics); % preallocates y points. Might need to change
colors = {[1,0,0],[0.9,0,0],[0.8,0,0],[0.7,0,0],[0.6,0,0],[0,1,0],[0,0.9,0],[0,0.8,0],[0,0.7,0],[0,0.6,0],[0,0,1],[0,0,0.9],[0,0,0.8],[0,0,0.7],[0,0,0.6],[1,0,1],[0.9,0,0.9],[0.8,0,0.8],[0.7,0,0.7],[0.6,0,0.6],[0,1,1],[0,0.9,0.9],[0,0.8,0.8],[0,0.7,0.7],[0,0.6,0.6],[1,1,0],[0.9,0.9,0],[0.8,0.8,0],[0.7,0.7,0],[0.6,0.6,0],[.85,0.325,0],[0.9,0.6,0.1],[0.4,0.2,0.6],[0.6,0.4,0.8],[0.3,0.74,0.9],[1,1,1],[0.8,0.8,0.8],[0.6,0.6,0.6]};

for i = 1:numrows
   figure(numpics+1); 
    for j = 1:numcols
    ypoints = {};
        for k = 1:numpics % need to change
        ypoints{k} = averages{i,j,k};      
        end
      if isnan(ypoints{k}) == 0 %determines if the data is good for graphing       
        ypoints = cell2mat(ypoints); 
        coefficients = polyfit(t,ypoints,3); % creates the coefficients of the fitted curve. degree changes
        newypoints = polyval(coefficients,t); % creates new y points that are smoooth
        if sqrt(numrows+1) - floor(sqrt(numrows+1))  <0.5
            p = 1;
        else 
            p = 0;
        end
        color = colors{j};
        subplot(ceil(sqrt(numrows+p)),ceil(sqrt(numrows)),i) % creates subplot
        plot(t,newypoints,'-','Color',color), hold on %can add real data points with t,y,'+'
        g =  plot(t,ypoints,'o','Color',color); %Adds true data points to graph
        gAnno = get(g, 'Annotation'); %following three lines make the true data points not appear in legend
        gLegend = get(gAnno, 'LegendInformation');
        set(gLegend,'IconDisplayStyle','off');
        title(['row ' num2str(i)])
        xlabel('time')
        ylabel('pixel value')
        ylim([ymin,ymax]); %specify y limits
        xlim([0,numpics-1]);
      end
    end
    if i == numrows
        for e = 1:35
            color = colors{e};
        f = color;
        g = color;
        subplot(ceil(sqrt(numrows+p)),ceil(sqrt(numrows)),i+1)
        h = plot(f,g,'Color',color);, hold on;
        legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35');       
        end
    end
end
%% Find the average slope of each square (ohgod)
for i = 1:numrows
    for j = 1:numcols
        for k = 1:numpics -1
            changesqaure{i,j,k}= averages{i,j,k+1} - averages{i,j,k};
        end
    end
end
avesquarechange = cell(numrows,numcols);
for i = 1:numrows
    for j =1:numcols
        avesquarechange{i,j} = nanmean(cell2mat(changesqaure(i,j,:)));
    end
end
%% Find the total change of each square
totsquarechange = cell(numrows,numcols);
for i = 1:numrows
    for j = 1:numcols
        totsquarechange{i,j} = averages{i,j,numpics} - averages{i,j,1}; % changed
%         if totsquarechange{i,j} >= 0
%             totsquarechange{i,j} = NaN;
%         end
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
% for i = 1:numrows
%     for j = 1:numcols
%         totchange{i,j} = (averages{i,j,numpics} - averages{i,j,1}); %creats a cell array of the change in T for each square
%         %if totchange{i,j} > 0
%          %   totchange{i,j} = NaN;
%         %end
%     end
% end

[xpoints,ypoints] = meshgrid(1:numcols,1:numrows);
zpoints = cell2mat(totsquarechange);
figure
surface(xpoints,ypoints,zpoints)
view(-37,64);
xlabel('Column')
ylabel('Row')
zlabel('Change in pixel value')
title('Total Temperature Change of sectios of Breast')
axis ij % makes axis match the figures
zpoints = cell2mat(avesquarechange);
figure
surface(xpoints,ypoints,zpoints)
view(-37,64);
xlabel('Column')
ylabel('Row')
zlabel('Average change in pixel value')
title('Average Temperature Change of sectios of Breast')
axis ij
 %% find the averages of each breast
answer = questdlg('Which side is the tumor on?','Tumor Side','Left','Right','Left')
figure(numpics)
prompt = ('Enter middle column');
dlgtitle = ('Middle column');
numlines = 1;

midcol = inputdlg(prompt,dlgtitle,numlines);

for i = 1:numrows % Takes data from right breast
    for j = 1:str2double(midcol{1})
        for k = 1:numpics % will have start at 1 
            Rbreastmean{i,j,k} = averages{i,j,k};
            if Rbreastmean{i,j,k} == 0
                Rbreastmean{i,j,k} = NaN;
            end
        end
    end 
end
totRbreastmean = cell(1,numpics ); % allocate arrays
Rstdv = cell(1,numpics);
for k =1:numpics % will have to get rid of -1. Finds averages of Rbreast and stdv 
    rightmean = cell2mat(Rbreastmean(:,:,k));
    rightmean = nanmean(rightmean);
    totRbreastmean{k} = nanmean(rightmean);
    Rightstdv = cell2mat(Rbreastmean(:,:,k));
    Rightstdv = nanstd(Rightstdv);
    Rstdv{k} = nanstd(Rightstdv);
end
for i = 1:numrows %gets the data for the left breast
    for j = str2double(midcol{1}):numcols
        for k = 1:numpics % will have start at 1 
            Lbreastmean{i,j,k} = averages{i,j,k};
            if Lbreastmean{i,j,k} == 0
                Lbreastmean{i,j,k} = NaN;
            end
        end
    end 
end
for i = 1:numrows %Makes all zero cells = NaN
    for j = 1:str2double(midcol{1})-1
        for k = 1:numpics % will have start at 1 
            Lbreastmean{i,j,k} = NaN;
            
        end
    end 
end
totLbreastmean = cell(1,numpics );
Lstdv = cell(1,numpics);
for k =1:numpics % will have to get rid of -1. calculates mean and stdv
    leftmean = cell2mat(Lbreastmean(:,:,k));
    leftmean = nanmean(leftmean);
    totLbreastmean{k} = nanmean(leftmean);
    Leftstdv = cell2mat(Lbreastmean(:,:,k));
    Leftstdv = nanstd(Leftstdv);
    Lstdv{k} = nanstd(Leftstdv);
end

 % get everything into a matrix
t = 0:numpics-1;
totRbreastmean = cell2mat(totRbreastmean);
totLbreastmean = cell2mat(totLbreastmean);
Lbreaststdv = cell2mat(Lstdv);
Rbreaststdv = cell2mat(Rstdv);
if isequal(answer,'Left') == 1
    tottumbreastmean(v,:) = totLbreastmean;
    tumstdv(v,:) = Lbreaststdv;
    normbreastmean(v,:) = totRbreastmean;
    normstdv(v,:) = Rbreaststdv;
elseif isequal(answer,'Right') == 1
    tottumbreastmean(v,:) = totRbreastmean;
    tumstdv(v,:) = Rbreaststdv;
    normbreastmean(v,:) = totLbreastmean;
    normstdv(v,:) = Lbreaststdv;
end

Rplot = errorbar(t,totRbreastmean,Rbreaststdv,'r');
hold on % plot everything
Lplot = errorbar(t,totLbreastmean,Lbreaststdv,'b');

title('Left vs Right Change Over Time')
legend('Right','Left')
xlabel('Time (min)')
ylabel('Pixel Value')
a = Rplot.Color;
Rplot.Color = 'r';
b = Lplot.Color;
Lplot.Color = 'b';

%% find average slope for each breast
for k = 1:numpics-1 % will have to make 1
    changeLbreast{k} =  totLbreastmean(k+1) - totLbreastmean(k);
end
aveLbreastchange = mean(cell2mat(changeLbreast));
for k = 1:numpics-1 % will have to make 1
    changeRbreast{k} =  totRbreastmean(k+1) - totRbreastmean(k);
end
aveRbreastchange = mean(cell2mat(changeRbreast));
%% Find total change for each breast
totLbreastchange = totLbreastmean(numpics) - totLbreastmean(1);
totRbreastchange = totRbreastmean(numpics) - totRbreastmean(1);
%% Graphs target area

% prompt = {'Top row','Bottom row','Left column','Right column'};
% dlgtitle = 'tumor region';
% answers = inputdlg(prompt,dlgtitle,1);
% warning('off')
% t = 0:numpics-1; % pictures start at t = 0, first pictue is lower than second. Strange
% ypoints = cell(1,numpics); % preallocates y points. Might need to change
% colors = {[1,0,0],[0.9,0,0],[0.8,0,0],[0.7,0,0],[0.6,0,0],[0,1,0],[0,0.9,0],[0,0.8,0],[0,0.7,0],[0,0.6,0],[0,0,1],[0,0,0.9],[0,0,0.8],[0,0,0.7],[0,0,0.6],[1,0,1],[0.9,0,0.9],[0.8,0,0.8],[0.7,0,0.7],[0.6,0,0.6],[0,1,1],[0,0.9,0.9],[0,0.8,0.8],[0,0.7,0.7],[0,0.6,0.6],[1,1,0],[0.9,0.9,0],[0.8,0.8,0],[0.7,0.7,0],[0.6,0.6,0],[.85,0.325,0],[0.9,0.6,0.1],[0.4,0.2,0.6],[0.6,0.4,0.8],[0.3,0.74,0.9]};
% toprow = str2double(answers{1}); 
% bottomrow = str2double(answers{2}); 
% leftcol = str2double(answers{3}); 
% rightcol = str2double(answers{4}); 
% for i = toprow:bottomrow
%    figure(numpics + 4); 
%     for j = leftcol:rightcol
%     ypoints = {};
%         for k = 1:numpics % need to change
%         ypoints{k} = averages{i,j,k};      
%         end
%       if isnan(ypoints{k}) == 0 %determines if the data is good for graphing       
%         ypoints = cell2mat(ypoints); legend('show')
%         coefficients = polyfit(t,ypoints,3); % creates the coefficients of the fitted curve. degree changes
%         newypoints = polyval(coefficients,t); % creates new y points that are smoooth
%         if sqrt(numrows) - floor(sqrt(numrows))  <0.5
%             k = 1;
%         else 
%             k = 0;
%         end
%         color = colors{j};
%         subplot(ceil(sqrt(bottomrow-toprow)),ceil(sqrt(bottomrow-toprow)),i-toprow+1) % creates subplot
%         r = rand;, g = rand;, b = rand; %sets random rgb values to make lots of colots
%         plot(t,newypoints,'-','Color',color,'DisplayName',num2str(j)), hold on %can add real data points with t,y,'+'
%         g =  plot(t,ypoints,'o','Color',color); %Adds true data points to graph
%         gAnno = get(g, 'Annotation'); %following three lines make the true data points not appear in legend
%         gLegend = get(gAnno, 'LegendInformation');
%         set(gLegend,'IconDisplayStyle','off');
%         title(['row ' num2str(i)])
%         xlabel('time')
%         ylabel('pixel value')
%         ylim([ymin,ymax]); %specify y limits
%         xlim([0,numpics-1]);
%       end
%     end
% end

%% Finds data for tumor and corresponding region
if answer == "Patient"
totsquarechange = cell2mat(totsquarechange);
[maximum, maxidx] = maxk(totsquarechange(:),10);
[lowrow, lowcol] = ind2sub(size(totsquarechange),maxidx);
warmregionidentifier = I_mat;
for k = 1:numpics
for i = 1:numel(lowcol)
    if lowrow(i) == 1 % fix error if first row or col = 1
        p = 1;
    elseif lowcol(i) == 1
        p = 1;
    else
        p = 0;
    end
   warmregionidentifier{k}(squareside*(lowrow(i)-1)+p:squareside-p:squareside*(lowrow(i)),...
       squareside*(lowcol(i)-1)+p:1:squareside*lowcol(i),:) = 10000;% converts every warmest regions to black
   warmregionidentifier{k}(squareside*(lowrow(i)-1)+p:1:squareside*(lowrow(i))...
       ,squareside*(lowcol(i)-1)+p:squareside-p:squareside*lowcol(i),:) = 10000; 
end
end
figure('Name','Crop Tumor Region')
[tumorcrop,tumrect] = imcrop(warmregionidentifier{numpics},[min(I_adj1) max(I_adj1)]); %sets the rectangle to crop all images
figure('Name','Crop corresponging region')
[corrcrop,corrrect] = imcrop(warmregionidentifier{numpics},[min(I_adj1) max(I_adj1)]); % sets crop for corresponding area
Itumor = cell(1,numpics);
Icorr = cell(1,numpics);
 for k =1:numpics % creates the tumor and corresponding area images
     Itumor{k} =imcrop(I_mat{k},tumrect);
     Icorr{k} = imcrop(I_mat{k}, corrrect); 
 end
 [tumrow,tumcol] = size(Itumor{1});
 for k = 1:numpics % removes black pixels from grid
    J = Itumor{k};
   J = double(J);
     for i = 1:tumrow
        for j = 1:tumcol
            if J(i,j) == 0 || J(i,j) == 10000
               J(i,j) = NaN;
            end
        end
     end
     Itumor{k} = J;
     
 end
 [corrrow,corrcol] = size(Icorr{1});
for k = 1:numpics
   K = Icorr{k};
   K = double(K);
    for i = 1:corrrow
        for j =1:corrcol
            if K(i,j) == 0 || K(i,j) == 10000
                K(i,j) = NaN;
            end
        end
    end
    Icorr{k} = K;
end
tumorregion = cell(1,numpics); %treats tumor and corresponding region as a single square
corrregion = cell(1,numpics);
for k = 1:numpics
    tumorregion{k} = nanmean(Itumor{k});
    corrregion{k} = nanmean(Icorr{k});
    tumorregion{k} = nanmean(tumorregion{k});
    corrregion{k} = nanmean(corrregion{k});
    tumstdv{k} = nanstd(Itumor{k});
    corrstdv{k} = nanstd(Icorr{k});
    tumstdv{k} = nanstd(tumstdv{k});
    corrstdv{k} = nanstd(corrstdv{k});
%     tumorregion{k} = nanmean(tumorregion{k});
%     corrregion{k} = nanmean(corrregion{k});
end
tumstdv = cell2mat(tumstdv);
corrstdv = cell2mat(corrstdv);

tottumorchange = tumorregion{numpics} - tumorregion{1}; %calculates the total change of the region
totcorrchange = corrregion{numpics} - corrregion{1};
figure
t = 0:numpics-1;
errorbar(t,totRbreastmean,Rbreaststdv,'r'), hold on
errorbar(t,totLbreastmean,Lbreaststdv,'b')
errorbar(t,cell2mat(tumorregion),tumstdv,'g')
errorbar(t,cell2mat(corrregion),corrstdv,'m')
legend('Right Breast','Left Breast','Tumor region','Corresponding region')
xlabel('Time (min)')
ylabel('Pixel Value')
title('Tumor region vs Corresponding Region Comapared to both Breasts')
changetumor = cell(1,numpics-1);
changecorr = cell(1,numpics-1);
for k = 1:numpics-1
    changetumor{k} = tumorregion{k+1} - tumorregion{k};
    changecorr{k} = corrregion{k+1} - corrregion{k};
end
avetumorchange = nanmean(cell2mat(changetumor));
avecorrchange = nanmean(cell2mat(changecorr));
figure
c = categorical({'Rbreast','Lbreast','Tumor region','Corresponding region'});
bardata = [totRbreastchange,totLbreastchange,tottumorchange,totcorrchange];
bar(c,bardata)
title('Total Change')
figure
bardata = [aveRbreastchange,aveLbreastchange,avetumorchange,avecorrchange];
bar(c,bardata)
title('Average Rate of change')
else
end
%% identiying regions of low change and highlighting them. Then comparing to ewach breast
if answer == "Volunteer"
totsquarechange = cell2mat(totsquarechange);
[maximum, maxidx] = maxk(totsquarechange(:),10);
[lowrow, lowcol] = ind2sub(size(totsquarechange),maxidx);
corrregionidentifier = I_mat;
for k = 1:numpics
for i = 1:numel(lowcol)
    if lowrow(i) == 1 % fix error if first row or col = 1
        p = 1;
    elseif lowcol(i) == 1
        p = 1;
    else
        p = 0;
    end
   corrregionidentifier{k}(squareside*(lowrow(i)-1)+p:squareside-p:squareside*(lowrow(i)),...
       squareside*(lowcol(i)-1)+p:1:squareside*lowcol(i),:) = 10000;% converts every warmest regions to black
   corrregionidentifier{k}(squareside*(lowrow(i)-1)+p:1:squareside*(lowrow(i))...
       ,squareside*(lowcol(i)-1)+p:squareside-p:squareside*lowcol(i),:) = 10000; 
end
end
figure('Name','Crop warm Region')
[warmcrop,warmrect] = imcrop(corrregionidentifier{numpics},[min(I_adj1) max(I_adj1)]); %sets the rectangle to crop all images
figure('Name','Crop Corresponding warm Region')
[corrwarmcrop,corrwarmrect] = imcrop(corrregionidentifier{numpics},[min(I_adj1) max(I_adj1)]); % sets crop for corresponding area
Iwarm = cell(1,numpics);
Icorrwarm = cell(1,numpics);
 for k =1:numpics % creates the tumor and corresponding area images
     Iwarm{k} =imcrop(corrregionidentifier{k},warmrect);
     Icorrwarm{k} = imcrop(corrregionidentifier{k}, corrwarmrect); 
 end
 [warmrow,warmcol] = size(Iwarm{1});
 for k = 1:numpics % removes white pixels from grid
    J = Iwarm{k};
   J = double(J);
     for i = 1:warmrow
        for j = 1:warmcol
            if J(i,j) == 10000 || J(i,j) == 0
               J(i,j) = NaN;
            end
        end
     end
     Iwarm{k} = J;
     
 end
 [corrwarmrow,corrwarmcol] = size(Icorrwarm{1});
for k = 1:numpics
   K = Icorrwarm{k};
   K = double(K);
    for i = 1:corrwarmrow
        for j =1:corrwarmcol
            if K(i,j) == 10000 || K(i,j) == 0
                K(i,j) = NaN;
            end
        end
    end
    Icorrwarm{k} = K;
end
warmregion = cell(1,numpics); %treats tumor and corresponding region as a single square
corrwarmregion = cell(1,numpics);
warmstdv = cell(1,numpics);
corrwarmstdv = cell(1,numpics);
for k = 1:numpics
    warmregion{k} = nanmean(Iwarm{k});
    corrwarmregion{k} = nanmean(Icorrwarm{k});
    warmregion{k} = nanmean(warmregion{k});
    corrwarmregion{k} = nanmean(corrwarmregion{k});
    warmstdv{k} = nanstd(Iwarm{k});
    corrwarmstdv{k} = nanstd(Icorrwarm{k});
    warmstdv{k} = nanstd(warmstdv{k});
    corrwarmstdv{k} = nanstd(corrwarmstdv{k});
%     tumorregion{k} = nanmean(tumorregion{k});
%     corrregion{k} = nanmean(corrregion{k});
end
warmstdv = cell2mat(warmstdv);
corrwarmstdv = cell2mat(corrwarmstdv);

totwarmchange = warmregion{numpics} - warmregion{1}; %calculates the total change of the region
totcorrwarmchange = corrwarmregion{numpics} - corrwarmregion{1};

t = 0:numpics-1;
errorbar(t,totRbreastmean,Rbreaststdv,'r'), hold on
errorbar(t,totLbreastmean,Lbreaststdv,'b')
errorbar(t,cell2mat(warmregion),warmstdv,'g')
errorbar(t,cell2mat(corrwarmregion),corrwarmstdv,'m')
legend('Right Breast','Left Breast','Warm region','Corresponding Warm region')
xlabel('Time (min)')
ylabel('Pixel Value')
title('Warm region vs Corresponding Warm Region Comapared to both Breasts')
changewarm = cell(1,numpics-1);
changecorrwarm = cell(1,numpics-1);
for k = 1:numpics-1
    changewarm{k} = warmregion{k+1} - warmregion{k};
    changecorrwarm{k} = corrwarmregion{k+1} - corrwarmregion{k};
end
avewarmchange = nanmean(cell2mat(changewarm));
avecorrwarmchange = nanmean(cell2mat(changecorrwarm));
figure
c = categorical({'Rbreast','Lbreast','Warm region','Corresponding Warm region'});
bardata = [totRbreastchange,totLbreastchange,totwarmchange,totcorrwarmchange];
bar(c,bardata)
title('Total Change')
figure
bardata = [aveRbreastchange,aveLbreastchange,avewarmchange,avecorrwarmchange];
bar(c,bardata)
title('Average Rate of change')
else
    
end

%% graphing total change across all patients
figure
t = 0:14;
alltumbreastmean = mean(tottumbreastmean);
alltumbreaststdv = mean(tumstdv);
allnormbreastmean = mean(normbreastmean);
allnormbreaststdv = mean(normstdv);

    tumplot = errorbar(t,alltumbreastmean,alltumbreaststdv,'r'); hold on
    normplot = errorbar(t,allnormbreastmean,allnormbreaststdv,'b');

title('Total tumor breast vs Normal breast')
legend('Tumor Breast','Normal Breast')
xlabel('Time (min)')
ylabel('Pixel Value')
a = tumplot.Color;
tumplot.Color = 'r';
b = normplot.Color;
normplot.Color = 'b';
