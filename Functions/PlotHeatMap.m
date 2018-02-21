function HMobj= PlotHeatMap(HeatMapData,RowLabels,ChannelsOut);
    granularity=21;
    colormap=redblue(granularity);
    try
        HMobj = clustergram(HeatMapData,'ColumnLabels',ChannelsOut,'RowLabels',RowLabels,'Standardize',1,'DisplayRange',3,'DisplayRatio',0.05,'Colormap',colormap,'ColumnPDist','correlation','ColumnLabelsRotate',45);
    catch
        close
        try
            HMobj = clustergram(HeatMapData,'ColumnLabels',ChannelsOut,'RowLabels',RowLabels,'Standardize',1,'DisplayRange',3,'DisplayRatio',0.05,'Colormap',colormap,'ColumnLabelsRotate',45);
        catch
            msgbox('Enter More Clusters to Make HeatMap');
        end
    end
end