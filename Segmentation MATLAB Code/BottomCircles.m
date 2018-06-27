[m,n] = size(circ_matrix);
c = 1;
d = 1;
flag = 0;
e = 0;

for i=round(m/2):m
    for j=round(n/2):n
        if circ_matrix(i,j) == 1
            if (circ_matrix(i-1,j) == 0)
                e = e+1;
                a = 1;
                b = 1;
                for k=i-2:i
                    for l=j-1:j+1
                        valuesAbove(a) = I(k,l);
                        a = a+1;
    %                     rectangle('Position',[l k 10 20])
                    end                              
                end
                meanA = mean(valuesAbove);            
                for k=i:i+2
                    for l=j-1:j+1
                        valuesBelow(b) = I(k,l);
                        b = b+1;
    %                     rectangle('Position',[l k 10 20])
                    end
                end
                meanB = mean(valuesBelow);
                difference = abs(meanA-meanB);
                dif(d) = difference; 
                d = d+1;
                if difference>100
                    points(c,1) = i;
                    points(c,2) = j;
                    c = c+1;
%                 else
%                      flag = 1;
%                      break                 
                 end
            end
        end
        if(flag==1)
          break
        end
    end
    if(flag==1)
      break
    end
end
A = [points(:,2);points(1,2)];
B = [points(:,1);points(1,1)];
figure; 
imshow(I,[])
line(Ap,Bp,'Color',[1 0 0],'linewidth',1) ;
getpts
[C,indC] = unique(points,'rows');
Ap = [C(:,2);C(1,2)];
Bp = [C(:,1);C(1,1)];
figure; 
imshow(I,[])
line(Ap,Bp,'Color',[1 0 0],'linewidth',1) ;
% plot(A,B)