function [HeatMapData,RowLabels,SizeCluster]=GetHeatMapData(num_clusters,idx,y,ClusterMethod,ClusterContrib);

     ClusterIter=[1:num_clusters];
    
    for j=ClusterIter;
            clusterselect=idx==j;
            SizeCluster(j)=sum(clusterselect);
            clusterselect2=y(clusterselect,:);
            if size(clusterselect2,1)==1
                HeatMapData(j,:)=clusterselect2;
            else
                HeatMapData(j,:)=median(clusterselect2); 
            end
            RowLabels{j}=strcat('Cluster ',num2str(j),' = ',num2str(ClusterContrib(j,3)),'%');
    end
end