function [newClustInfo, maxNumClusters] = symmetry_alignClusters(ClustInfo)
newClustInfo = ClustInfo; %Copy CLUSTINFO 
% ClustCentroids = {};
CentInd = [];

for a = 1:15
   ImStruct = ClustInfo{a,1};
   for b = 1:length(ImStruct)
       if (isempty(CentInd))
           CentInd(1,:) = ImStruct(b).ClusterCentroid;
       else
           CentInd(end+1,:) = ImStruct(b).ClusterCentroid;
       end
   end
end
% DBSCAN ON CENTROIDS
epsilon = 15; %Best Epsilon by Testing, maybe change per patient? GUI Option
minPts = 0;

CentIDX = DBSCAN(CentInd,epsilon,minPts);

spt = 1;
for a = 1:15
    for b = 1:length(newClustInfo{a,1})
        newClustInfo{a,1}(b).NormalizedCluster = CentIDX(spt);
        spt=spt+1;
    end
end

maxNumClusters = max(CentIDX);

end