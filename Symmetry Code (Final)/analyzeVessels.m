%Look at vessels

clearvars -except I_mat xR yR xL yL,
close all, 
%% Location input and Nipple selection
%Only redo if the patient/volunteer has changed
redo = questdlg('New Patient/Volunteer?/Different Side', 'Cancel', 'No');
if strcmp(redo, 'Yes') == 1%% Input
    path  = uigetdir;
    location = strcat(path, '\');
    a=0;
    for i=1:15          
        I_mat{i} = imread([location sprintf('%04d.tif',a)]);    % Read each image into I_mat
        a=a+120;            % Go to next image (for cropped)
    end

%     rightOrLeft = questdlg('Right Side?', 'Yes', 'No');
%     if strcmp(rightOrLeft, 'Yes') == 1
%         figure('Name','Select nipple (Right)'), 
%          for i = 1:15
%             I1 = I_mat{i}(find(I_mat{i}>0));
%             imshow(I_mat{i}, [min(I1) max(I1)])
%             hold on
%             [xR{i},yR{i}] = ginput(1);
%          end
%          for i = 1:15
%             xMovement{i} = xR{i} - xR{7};
%             yMovement{i} = yR{i} - yR{7};
%         end
%     else
%         figure('Name','Select nipple (Left)'), 
%         for i =1:15
%             I1 = I_mat{i}(find(I_mat{i}>0));
%             imshow(I_mat{i}, [min(I1) max(I1)]) % gets coordinates of nipple
%             hold on,
%             [xL{i},yL{i}] = ginput(1);
%         end
%         for i = 1:15
%             xMovement{i} = xL{i} - xL{7};
%             yMovement{i} = yL{i} - yL{7};
%         end
%     end
end

for i = 1:15
    figure,
    I1 = I_mat{i}(find(I_mat{i}>0));
    imshow(I_mat{i}, [min(I1) max(I1)])
    rect1 = imrect();
    position{i} = wait(rect1);
end 


for i = 1:15
    I_cropped{i} = imcrop(I_mat{i}, position{i});
    avg(i) = mean2(I_cropped{i});
end 

for i = 1:14
    stepChange(i) = avg(i+1)- avg(i);
end 

avgStepChange = mean2(stepChange);
totalChange = avg(15) - avg(1);


