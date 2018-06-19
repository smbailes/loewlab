function [location, ptID] = choosedialog
%% Patient or Volunteer
    answer = questdlg('ID Patient Type:','Patient Type','Patient','Volunteer','Patient');
    if (strcmp(answer, 'Patient'))
        ptID = patientselect; %return patient ID
    else
        ptID = volunteerselect; %return patient ID
    end
    
%% File Explorer 
    dir = uigetdir; %file exploerer to find image directory
    c = questdlg('Use Cropped Images?', 'Yes', 'No'); %check to see if you want to use cropped images
    if(strcmp(c, 'Yes'))
        location = strcat(dir, '\Cropped\');
    else
        location = strcat(dir, '\');
    end
end