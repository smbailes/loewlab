clear; 
close all;
for k=1:2
    I = imread('1799 - V10.tif');
    I = im2double(I); 
    figure, imshow(I,[]); [y,x] = getpts;
    P = [x(:) y(:)];
    Options = struct;
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

    [O,J] = Snake2D(I,P,Options);
%     figure;
%     imshow(I,[])
%     hold on
%     plot([O(:,2);O(1,2)],[O(:,1);O(1,1)]);
    
    Y = immultiply(J,I);

    [m,n] = size(I);
    for i=1:m
        for j=1:n
            if Y(i,j) == 0
                Y(i,j) = min(min(I));                
            end
            A(i,j,k) = Y(i,j);
        end
    end      
%     figure;
%     imshow(Y,[])
    O2(:,:,k) = O(:,:); 
end
figure;
imshow(I,[])
hold on
plot([O2(:,2,1);O2(1,2,1)],[O2(:,1,1);O2(1,1,1)]);
hold on
plot([O2(:,2,2);O2(1,2,2)],[O2(:,1,2);O2(1,1,2)]);
H = imadd(A(:,:,1),A(:,:,2));
figure; imshow(H,[])