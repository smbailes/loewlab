function I = imread_irfanview(filename)
% This function IMREAD_IRFANVIEW is to load images from disk which
% have a format not supported by Matlab, but support by Irfanview
%
% If irfanview is not found, adjust the path inside the function
%

% Irfanview location - Put the write directory based on the computer you're
% using

%cm = '"C:\Users\kamona\Documents\IrfanView\i_view32.exe"'; %LAUTERBUR
%cm = '"C:\Program Files (x86)\IrfanView\i_view32.exe"'; %LOVELACE
%cm = '"C:\Program Files\IrfanView\i_view64.exe"'; %laptop directory
cm = ' "C:\Users\Zainab Mahmood\Documents\MATLAB\IrfanView\i_view32.exe" ';


% Make temporary filename
temp_folder = getenv('temp');
filename_temp =[temp_folder  '\test'  sprintf('%08.0f',now*100) '.tif'];

% Command line command
cm = [cm ' ' filename '  /convert=' filename_temp];

% Execute command line command
[status,result] = system(cm);

if(status~=0) 
    error(result)
    I = [];
else 
    I =imread(filename_temp);
    delete(filename_temp);
end