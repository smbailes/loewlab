close all;

edgeprewitt = edge(I,'prewitt');
figure,imshow(edgeprewitt)
title('Prewitt edges');

edgesobel = edge(I,'sobel');
figure,imshow(edgesobel)
title('Sobel edges');

edgeroberts = edge(I,'roberts');
figure,imshow(edgeroberts)
title('Roberts Edges')

edgeslog = edge(I,'log');
figure,imshow(edgeslog)
title('Log Edges');