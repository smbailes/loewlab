function userselect = choosedialog
    d = dialog('Position',[300 300 250 150],'Name','Select Patient');
    
    txt = uicontrol('Parent',d,'Style','text','Position',[20 80 210 40],... 
           'String','Select a User');
       
    popup = uicontrol('Parent',d,'Style','popup','Position',[75 70 100 25],...
           'String',{'Sydney';'Samhita';},... //ADD NEW USERS HERE
           'Callback',@popup_callback);
       
    btn = uicontrol('Parent',d,'Position',[89 20 70 25],... % button
           'String','Next','Callback','delete(gcf)');     
    
    userselect = 'Sydney';   %Default Answer, should be set to first name on list
    
    % Wait for d to close before running to completion
    uiwait(d);
       function popup_callback(popup,event)
          idx = popup.Value;
          popup_items = popup.String;
          userselect = char(popup_items(idx,:));
       end
end