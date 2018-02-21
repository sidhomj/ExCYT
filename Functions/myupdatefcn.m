function txt = myupdatefcn(empt,event_obj,idx,Y)
    pos = get(event_obj,'Position');
    loc=find(Y==pos);
    cluster=idx(loc(1));
    txt=strcat('Cluster_', num2str(cluster));
end