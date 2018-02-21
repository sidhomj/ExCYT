function [num,ChannelsOut]= ReadCSVFile(fileread1)
    
    %fid=fopen(fileread1);
    [fid,msg]=fopen(fileread1);
    
    if fid<0
        disp(msg)
    end
    
    header=textscan(fid,repmat('%s',1,100),1,'delimiter',',');
    data=textscan(fid,repmat('%f',1,100),'delimiter',',');

    for i=1:size(data,2);
        try
       data2(:,i)=data{1,i};
        catch
            continue
        end
    end

    try
        channelcheck=data2(1,:);
        channelcheck=isnan(channelcheck);
        channelcheck=find(channelcheck==1);
        Num_Channels=channelcheck(1)-1;
    catch
        Num_Channels=size(data2,2);
    end

    Channels=header(1:Num_Channels);
    num=data2(:,[1:Num_Channels]);

    for i=1:size(Channels,2);
        temp=Channels{i};
        temp=temp{1};
        %temp=strsplit(temp,':: ');
        ChannelsOut{i}=temp;
    end
    
end
