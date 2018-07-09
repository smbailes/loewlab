% clear; 
% close all;

% I = imread('1799 - V10.tif');
% I = im2double(I); 
figure, imshow(I,[]); [y,x] = getpts;
[m,n] = size(total);
%y = [0; 0.5*m; m; m; m; m; m; 0.5*m; 0];
%x = [n; n; n; 0.75*n; 0.5*n; 0.25*n; 0; 0; 0];

% k = 0; m = 0;
% for i = 1:1:length(x)
%     k = k + 1;
%     xx(k,1) = x(i,1);
%     yy(k,1) = y(i,1);
% end

P = [x(:),y(:)];

figure;
imshow(I,[])
hold on
set(displ, 'AlphaData', total)
hold off;
Options = struct;

%% Basic
Options.Verbose = true;
Options.Iterations = 300;
Options.nPoints = 300;
[O,J] = Snake2D(connectedtop,P,Options);

%% GVF
clc
Options.Verbose = true;
Options.nPoints = 500;
Options.Gamma = 2;
Options.Iterations = 400;

Options.Sigma1 = 30;
Options.Wline = 0.01;
Options.Wedge = 25;
Options.Wterm = 0.01;
Options.Sigma2 = 5;

Options.Mu = 0.02;
Options.GIterations = 600;
Options.Sigma3 = 3;

Options.Alpha = 1.0;
Options.Beta = 0.5;
Options.Delta = 0.00050;
Options.Kappa = 4;

[O,J] = Snake2D(connectedtop,P,Options);
sprintf('done')
  %% Result
figure;
imshow(I,[])
hold on
plot([O(:,2);O(1,2)],[O(:,1);O(1,1)]);

