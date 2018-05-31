function [K]=circ(grote)

K=zeros(300,300)+100;

if grote>10

grote=10;

end

for radius=0:.05:grote

theta=0:.01:7;

[X,Y] = pol2cart(theta,radius);

A=[round(X)+11]; B=[round(Y)+11];

ring = sub2ind(size(K),A,B); K(ring) = 190;

end

imshow(K, [0 255]);

end