% clear; 
% close all;

% I = imread('1799 - V10.tif');
% I = im2double(I); 
%figure, imshow(I,[]); [y,x] = getpts;
[m,n] = size(total);
y = [0.25*m; 0.5*m; 0.8*m; 0.9*m; 0.7*m; 0.9*m; 0.8*m; 0.5*m; 0.25*m];
x = [n; n; n; 0.75*n; 0.5*n; 0.25*n; 0; 0; 0];
% [y,x] = find(total == 65536);
% 
% k = 0;
% for i = length(x):-1:1
%     k = k + 1;
%     xx(k,1) = x(i,1);
%     yy(k,1) = y(i,1);
% end

P = [y(:),x(:)];

figure;
imshow(I,[])
hold on
plot(P)
%set(displ, 'AlphaData', total)
hold off;
title('Image for Snakes');
Options = struct;

%% Basic
Options.Verbose = true;
Options.Iterations = 300;
Options.nPoints = 300;
[O,J] = Snake2D(I,P,Options);

%% GVF
clc
Options.Verbose = true;
Options.nPoints = 600;
Options.Gamma = 2;
Options.Iterations = 3000;

Options.Sigma1 = 20;
Options.Wline = 0.06;
Options.Wedge = 25;
Options.Wterm = 0.01;
Options.Sigma2 = 5;

Options.Mu = 0.02;
Options.GIterations = 600;
Options.Sigma3 = 3;

Options.Alpha = 2;
Options.Beta = 1;
Options.Delta = 0.0050;
Options.Kappa = 4;

[O,J] = Snake2D(total,P,Options);
sprintf('done')
  %% Result
figure;
imshow(I,[])
hold on
plot([O(:,2);O(1,2)],[O(:,1);O(1,1)]);

