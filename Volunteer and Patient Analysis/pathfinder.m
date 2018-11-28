function [location, ptID,answer] = choosedialog
%% Patient or Volunteer
    answer = questdlg('ID Patient Type:','Patient Type','Patient','Volunteer','Patient');
    if (strcmp(answer, 'Patient'))
        ptID = patientselect; %return patient ID
    else
        ptID = volunteerselect; %return patient ID
    end
    
%% File Explorer 
    dir = uigetdir; %file exploerer to find image directory

    location = strcat(dir, '\');

end