function [G,GraphStore]=CreateGraph(data,k);

clear G
distance = 'euclidean';
[IDX,D] = knnsearch(data,data,'k',k+1,'distance',distance);
IDX(:,1) = []; D(:,1) = [];

G = knn2jaccard(IDX);
G=G+G';
GraphStore=graph(G);

end