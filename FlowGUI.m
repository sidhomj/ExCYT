function varargout = FlowGUI(varargin)
% FLOWGUI MATLAB code for FlowGUI.fig
%      FLOWGUI, by itself, creates a new FLOWGUI or raises the existing
%      singleton*.
%
%      H = FLOWGUI returns the handle to a new FLOWGUI or the handle to
%      the existing singleton*.
%
%      FLOWGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLOWGUI.M with the given input arguments.
%
%      FLOWGUI('Property','Value',...) creates a new FLOWGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FlowGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FlowGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FlowGUI

% Last Modified by GUIDE v2.5 26-Sep-2017 17:16:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FlowGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FlowGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FlowGUI is made visible.
function FlowGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FlowGUI (see VARARGIN)
warning('off','all');
% Choose default command line output for FlowGUI
handles.output = hObject;
handles.clustermethodbox.String={'Hard KMEANS (on t-SNE)','Hard KMEANS (on HD Data)','DBSCAN','Hierarchical Clustering','Network Graph-Based','Self Organized Map','GMM - Expectation Minimization','Variational Bayesian Inference for GMM'};
addpath('Functions/');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FlowGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FlowGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fileread1_Callback(hObject, eventdata, handles)
% hObject    handle to fileread1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileread1 as text
%        str2double(get(hObject,'String')) returns contents of fileread1 as a double
if strcmp(handles.fileread1.String(end-3:end),'.csv')
    [num,ChannelsOut]= ReadCSVFile(handles.fileread1.String);
    handles.num=num;
    handles.ChannelsAll=ChannelsOut;
    handles.channelselect.String=ChannelsOut;
    handles.xaxis.String=ChannelsOut;
    handles.yaxis.String=ChannelsOut;
elseif strcmp(lower(handles.fileread1.String(end-3:end)),'.fcs')
    %[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(handles.fileread1.String);
    try
        [data, marker_names, channel_names, scaled_data, compensated_data, fcshdr] = readfcs_v2(handles.fileread1.String);
        data=transpose(data);
        compensated_data=transpose(compensated_data);
        handles.num=compensated_data;
        handles.num2=data;
        header=marker_names;       
    catch
        [data, fcshdr, fcsdatscaled, compensated_data] = fca_readfcs(handles.fileread1.String);
        data=double(data);
        headerdata=fcshdr.par;
        for i=1:size(headerdata,2);
            if ~isempty(headerdata(i).name2)
                header{i}=headerdata(i).name2;
            else
                header{i}=headerdata(i).name;
            end
        end
        handles.num=data;
        handles.num2=data;
    end
%     mindata=min(fcsdat);
%     for i=1:size(mindata,2);
%         if mindata(i)<0
%             mindataapp=repmat(mindata(i),size(fcsdat,1),1);
%             fcsdat(:,i)=fcsdat(:,i)+abs(mindataapp);
%         end
%     end
%     mindata=repmat(mindata,size(fcsdat,1),1);
%     fcsdat=fcsdat+abs(mindata);

    handles.ChannelsAll=header;
    handles.channelselect.String=header;
    handles.xaxis.String=header;
    handles.yaxis.String=header;
    handles.heatmaptsne.String=header;
end
guidata(hObject,handles);
msgbox('Data Imported');



% --- Executes during object creation, after setting all properties.
function fileread1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileread1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectFileButton.
function SelectFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=ResetGUI(handles);

[filename,folder]=uigetfile({'*.fcs'},'Select file','MultiSelect','on');
handles.files.String=filename;

if ~iscell(filename)
    filename=cellstr(filename);
end

if strcmp(handles.eventperfile.String,'Events Per File') || isempty(handles.eventperfile.String)
    subsample=0;
else
    subsample=1;
    sampleamnt=str2num(handles.eventperfile.String);
end

num=[];
num2=[];
for i=1:size(filename,2)
    filenamequery=filename(i);
    filenamequery=fullfile(folder,filenamequery);
    try
        [data, marker_names, channel_names, scaled_data, compensated_data, fcshdr] = readfcs_v2(filenamequery{1});
        data=single(transpose(data+1));
        compensated_data=single(transpose(compensated_data+1));
        if subsample==1
            if size(data,1)>sampleamnt || size(compensated_data,1)>sampleamnt
                data=datasample(data,sampleamnt);
                compensated_data=datasample(compensated_data,sampleamnt);
            end    
            %data=asinh(data/5);
            %compensated_data=asinh(data/5);
        end
        num=[num;compensated_data];
        num2=[num2;data;];
%         handles.num=compensated_data;
%         handles.num2=data;
        header=marker_names;       
    catch
        [data, fcshdr, fcsdatscaled, compensated_data] = fca_readfcs(filenamequery{1});
        data=single(double(data+1));
        if subsample==1
            if size(data,1)>sampleamnt
                data=datasample(data,sampleamnt);
            end 
             %data=asinh(data/5);
        end
        headerdata=fcshdr.par;
        for j=1:size(headerdata,2);
            if ~isempty(headerdata(j).name2)
                header{j}=headerdata(j).name2;
            else
                header{j}=headerdata(j).name;
            end
        end
%         handles.num=data;
%         handles.num2=data;
        num=[num;data];
        num2=[num2;data];
    end
end

handles.num=num; 
handles.num2=num2;
handles.ChannelsAll=header;
handles.channelselect.String=header;
handles.xaxis.String=header;
handles.yaxis.String=header;
handles.heatmaptsne.String=header;

guidata(hObject,handles);
msgbox('Data Imported');



%[filename,folder]=uigetfile({'*.fcs'},'Select file');
%filename=fullfile(folder, filename);
%set(handles.fileread1,'string',filename);
%fileread1_Callback(hObject,eventdata,handles);

    function handles=ResetGUI(handles);
        fieldremove={'num','num2','ChannelsAll','num_samples','channel_select','y2','ChannelsOut','y','Y','tsne_xlim','tsne_ylim'};
        for i=1:size(fieldremove,2);
            if isfield(handles,fieldremove{i})
                handles=rmfield(handles,fieldremove{i});
            end
        end
        
        handles.numsamples.String='';
        handles.perc_file.String='';
        handles.channelselect.Value=[];
        handles.channelselect.String={};
        handles.heatmaptsne.Value=1;
        handles.heatmaptsne.String={'Channel'};
        handles.popupmenu1.Value=1;
        handles.popupmenu1.String={'Channel'};
        handles.popupmenu2.Value=1;
        handles.popupmenu2.String={'Channel'};
        
        
        


% --- Executes on button press in tnsebutton.
function tnsebutton_Callback(hObject, eventdata, handles)
% hObject    handle to tnsebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num=handles.num;
if ~isfield(handles,'num_samples')
    msgbox('Number of Samples Not Entered','Error','error');
end

if ~isfield(handles,'y')
    y=datasample(num,handles.num_samples);
    handles.y2=y;
    if handles.inst_type.Value==1
        y=asinh(y/150);
        handles.transy2=y;
    elseif handles.inst_type.Value==2;
        y=asinh(y/5);
        handles.transy2=y;
    end
    
else
    if handles.num_samples>size(handles.y2,1)
        y=datasample(num,handles.num_samples);
        handles.y2=y;
        if handles.inst_type.Value==1
            y=asinh(y/150);
            handles.transy2=y;
        elseif handles.inst_type.Value==2;
            y=asinh(y/5);
            handles.transy2=y;
        end
    else
        y=handles.y2;
        y=datasample(y,handles.num_samples);
        handles.y2=y;
         if handles.inst_type.Value==1
            y=asinh(y/150);
            handles.transy2=y;
        elseif handles.inst_type.Value==2;
            y=asinh(y/5);
            handles.transy2=y;
        end
    end
end

if ~isfield(handles,'channel_select')
    msgbox('Select Channels For Analysis','Error','error');
end

y=y(:,handles.channel_select);
handles.y=handles.y2(:,handles.channel_select);
handles.transy=y;

handles.ChannelsOut=handles.ChannelsAll(handles.channel_select);
handles.popupmenu1.String=handles.ChannelsAll;
handles.popupmenu2.String=handles.ChannelsAll;

normalizetsne=1;
hbox=msgbox('Running t-SNE Analysis');
Y=tsne(y,'Standardize',normalizetsne);
close(hbox);
handles.Y=Y;

scatter(handles.axes1,Y(:,1),Y(:,2),'filled');
handles.axes1.XTickLabel={};
handles.axes1.YTickLabel={};
handles.axes1.Title.String='TSNE';
handles.tsne_xlim=handles.axes1.XLim;
handles.tsne_ylim=handles.axes1.YLim;

guidata(hObject,handles);



function clusterparameter_Callback(hObject, eventdata, handles)
% hObject    handle to clusterparameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clusterparameter as text
%        str2double(get(hObject,'String')) returns contents of clusterparameter as a double

key = get(gcf,'CurrentKey');
if(strcmp (key , 'return'))
    clusterbutton_Callback(hObject, eventdata, handles)
end



% --- Executes during object creation, after setting all properties.
function clusterparameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterparameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in clustermethodbox.
function clustermethodbox_Callback(hObject, eventdata, handles)
% hObject    handle to clustermethodbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clustermethodbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clustermethodbox

sel=handles.clustermethodbox.Value;
if ismember(sel,[1 2 6 7 8])
    set(handles.clusterparameter,'String','# of Clusters');
elseif ismember(sel,[3 4])
    set(handles.clusterparameter,'String','Distance Factor');
elseif ismember(sel,[5])
    set(handles.clusterparameter,'String','k-nearest neighbors');
end
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function clustermethodbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clustermethodbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numsamples_Callback(hObject, eventdata, handles)
% hObject    handle to numsamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numsamples as text
%        str2double(get(hObject,'String')) returns contents of numsamples as a double

handles.num_samples=str2num(handles.numsamples.String);
handles.perc_file.String=num2str(100*handles.num_samples/size(handles.num,1));
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function numsamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numsamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function perc_file_Callback(hObject, eventdata, handles)
% hObject    handle to perc_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of perc_file as text
%        str2double(get(hObject,'String')) returns contents of perc_file as a double

handles.num_samples=round((str2num(handles.perc_file.String)/100)*size(handles.num,1));
handles.numsamples.String=num2str(handles.num_samples);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function perc_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to perc_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in channelselect.
function channelselect_Callback(hObject, eventdata, handles)
% hObject    handle to channelselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channelselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channelselect
handles.channel_select=handles.channelselect.Value;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function channelselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channelselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in clusterbutton.
function clusterbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clusterbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fieldremove={'ManualClusterCount','colorspec1','idx','ClusterMethod','ClusterContrib','num_clusters','HeatMapData','RowLabels','SizeCluster','I','Imod','Ifinal','line','thresholdcurrent','thresholdbook','threshold_count','graphclustermethod'};

for i=1:size(fieldremove,2);
    if isfield(handles,fieldremove{i})
        handles=rmfield(handles,fieldremove{i});
    end
end
handles.threshlist.String={};
handles.cluster_sel.String={};
handles.cluster_sel.Value=[];
handles.clusterplot.String={};
handles.clusterplot.Value=[];
handles.clusterfreq.String='';
handles.popupmenu1.Value=1;
handles.popupmenu2.Value=1;

Y=handles.Y;
ClusterMethod=handles.clustermethodbox.Value;
clusterparameter=str2num(handles.clusterparameter.String);

switch ClusterMethod

    case 1
        
        hbox=msgbox('Clustering Events...');
        num_clusters=clusterparameter;
        idx=kmeans(Y,num_clusters,'Start','uniform');
    case 2
        
        hbox=msgbox('Clustering Events...');
        num_clusters=clusterparameter;
        idx=kmeans(handles.transy,num_clusters,'Start','uniform');
%     elseif ClusterMethod==3
%         
%         hbox=msgbox('Clustering Events...');
%         [centers,U] = fcm(Y,clusterparameter);
%         thresh=0.75;
%         parfor i=1:size(U,2);
%             probdist=U(:,i);
%             probdist=probdist>thresh
%             idxi=find(probdist);
%             if sum(idxi)==0
%                 idx(i)=clusterparameter+1;
%             else
%                 idx(i)=idxi;
%             end
%         end
%         num_clusters=clusterparameter+1;
%         idx=transpose(idx);
%     elseif ClusterMethod==4
%         hbox=msgbox('Clustering Events...');
%         [centers,U] = fcm(handles.y,clusterparameter);
%         thresh=0.75;
%         parfor i=1:size(U,2);
%             probdist=U(:,i);
%             probdist=probdist>thresh
%             idxi=find(probdist);
%             if sum(idxi)==0
%                 idx(i)=clusterparameter+1;
%             else
%                 idx(i)=idxi;
%             end
%         end
%         num_clusters=clusterparameter+1;
%         idx=transpose(idx);
        
    case 3
        hbox=msgbox('Clustering Events...');
        epsilonf=clusterparameter/100;
        D=pdist(Y);
        epsilon=(epsilonf)*median(D); %.02 default
        MinPoints=1;%(0.0001)*size(Y,1); %.0001 default
        [idx,isnoise]=DBSCAN(Y,epsilon,MinPoints);
        num_clusters=max(idx);
    case 4
        hbox=msgbox('Clustering Events...');
        dm=pdist(handles.transy);
        z=linkage(dm);
        idx=cluster(z,'cutoff',clusterparameter);
        num_clusters=max(idx);
    case 5
        NetworkGui(handles);
        waitfor(findobj('Tag','networkgui'));
        
        hbox=msgbox('Creating Graph...');
        [G,GGraph]=CreateGraph(handles.transy,clusterparameter);
        close(hbox);
        
        handles=guidata(findobj('Tag','clusterbutton'));
     
        switch handles.graphclustermethod  
            case 1      
                hbox=msgbox('Clustering Events...');
                N=length(G);
                W=PermMat(N);                     % permute the graph node labels
                A=W*G*W';
                
                %%RG Version
%                 [row,col,val]=find(A);
%                 input=table(row,col,val);
%                 writetable(input,'G.csv');
%                 RG_clust('G.csv')
                
                [COMTY ending] = cluster_jl_cppJW(A,1);
                J=size(COMTY.COM,2);
                VV=COMTY.COM{J}';
                idx=W'*VV;      
                              
            case 2
                hbox=msgbox('Clustering Events...');
                idx=GCModulMax2(G); 
            case 3
                hbox=msgbox('Clustering Events...');
                idx=GCModulMax3(G);
            case 4
                hbox=msgbox('Clustering Events...');
                idx=GCDanon(G);
            case 5
                clusterparameter2=inputdlg('Enter # of Clusters');
                hbox=msgbox('Clustering Events...');
                clusterparameter2=str2num(clusterparameter2{1});
                idx=GCSpectralClust1(G,clusterparameter2);
                idx=idx(:,clusterparameter2);
                
        end
        num_clusters=max(idx);
    case 6
        hbox=msgbox('Clustering Events...');
        net=selforgmap([round(sqrt(clusterparameter)),round(sqrt(clusterparameter))]);
        net.trainParam.showWindow = false;
        net=train(net,transpose(handles.transy));
        idx=transpose(vec2ind(net(transpose(handles.transy))));
        num_clusters=max(idx);
    case 7
        hbox=msgbox('Clustering Events...');
        try 
         idx=transpose(mixGaussEm(transpose(handles.transy),clusterparameter));
        num_clusters=max(idx); 
        catch
            msgbox('Enter smaller # of Clusters');
        end
        
    case 8
        hbox=msgbox('Clustering Events...');
        try
        idx=transpose(mixGaussVb(transpose(handles.transy),clusterparameter));
        num_clusters=max(idx);
        catch
            msgbox('Enter smaller # of Clusters');
        end
        
        end
    close(hbox);
    
    [colorspec1,colorspec2]=CreateColorTemplate(num_clusters);
    handles.colorspec1=colorspec1;
    
    clear colorscheme
    for i=1:size(idx,1);
        colorscheme(i,:)=colorspec1(idx(i)).spec;
    end
    scatter(handles.axes1,Y(:,1),Y(:,2),[],colorscheme,'filled');
    handles.axes1.XTickLabel={};
    handles.axes1.YTickLabel={};
    
    
    y=handles.transy;
    ClusterContrib=tabulate(idx);
    
    for i=1:num_clusters
        ClusterNames(i)=strcat('Cluster ',num2str(i),{' - '},num2str(ClusterContrib(i,3)),{'%'});
    end

    [HeatMapData,RowLabels,SizeCluster]=GetHeatMapData(num_clusters,idx,y,ClusterMethod,ClusterContrib);
    
    handles.cluster_sel.String=ClusterNames;
    handles.idx=idx; 
    handles.ClusterContrib=ClusterContrib;
    handles.num_clusters=num_clusters;
    handles.HeatMapData=HeatMapData;
    handles.RowLabels=RowLabels;
    handles.SizeCluster=SizeCluster;
    handles.I=[1:num_clusters];
    handles.Imod=handles.I;
    guidata(hObject,handles);
    datacursormode on
    dcm_obj=datacursormode(handles.axes1.Parent);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,idx,Y})
    
    function NetworkGui(handles)
        
         fh = figure('units','pixels',...
                  'units','normalized',...
                  'position',[0.25 0.25 .175 .125],...
                  'menubar','none',...
                  'name','Select Graph Clustering Method',...
                  'tag','networkgui',...
                  'numbertitle','off',...
                  'resize','off');
        guidata(fh,handles);
    
              
        clusteroptions=uicontrol('Style','listbox',...
        'String',{'Modularity Max - Louvain';'Modularity Max - Fast Greedy';'Modularity Max - Newman';...
        'Danon Method';'Spectral Clustering'},...
        'units','normalized',...
        'Position',[.1,.3,0.8,0.6],...
        'Tag','graphclusteropt');
    
        selectbutton=uicontrol('Style','pushbutton',...
            'String','Select',...
            'units','normalized',...
            'Position',[0.75,.1,.2,0.2],...
            'Tag','selectgraphcluster',...
            'Callback',@selectgraphcluster);
        
        function selectgraphcluster(hObject, eventdata)
            temp=findobj('Tag','graphclusteropt');
            clustergraphmethod=temp.Value;
            handles=guidata(hObject);
            handles.graphclustermethod=clustergraphmethod;
            guidata(findobj('Tag','clusterbutton'),handles);
            closereq;

     
    


% --- Executes on button press in selcluster.
function selcluster_Callback(hObject, eventdata, handles)
% hObject    handle to selcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dcm_obj=datacursormode(handles.axes1.Parent);
set(dcm_obj,'Enable','off');

if ~isfield(handles,'ManualClusterCount')
    if isfield(handles,'idx')
        handles=rmfield(handles,'idx');
        handles=rmfield(handles,'colorspec1');
    end
    [colorspec1,colorspec2]=CreateColorTemplate(100);
    handles.colorspec1=colorspec1;
    handles.cluster_sel.String={};
    Y=handles.Y;
    scatter(handles.axes1,Y(:,1),Y(:,2),'filled');
    handles.axes1.XTickLabel={};
    handles.axes1.YTickLabel={};
    handles.axes1.Title.String='TSNE';
    hold(handles.axes1)
    handles.ManualClusterCount=1;
    I=1;
else
    hold(handles.axes1)
    handles.ManualClusterCount=handles.ManualClusterCount+1;
    colorspec1=handles.colorspec1;
    idx=handles.idx;
    I=handles.I;
    I=[I,I(end)+1];
    handles.I=I;
    handles.Imod=I;
    ClusterContrib=handles.ClusterContrib;
end


sel=selectdata('SelectionMode','Lasso','Verify','on');
if handles.ManualClusterCount~=1
    sel=sel{handles.ManualClusterCount};
end
in=zeros(size(handles.Y,1),1);
in(sel)=1;
in=logical(in);
%[x,y]=ginput;
%in=inpolygon(handles.Y(:,1),handles.Y(:,2),x,y);
Yplot=handles.Y(in,:);
hold(handles.axes1);
scatter(handles.axes1,Yplot(:,1),Yplot(:,2),[],colorspec1(handles.ManualClusterCount).spec,'filled');
handles.axes1.XLim=handles.tsne_xlim;
handles.axes2.YLim=handles.tsne_ylim;
hold(handles.axes1);

ClusterContrib(handles.ManualClusterCount,1)=handles.ManualClusterCount;
ClusterContrib(handles.ManualClusterCount,2)=sum(in);
ClusterContrib(handles.ManualClusterCount,3)=100*(sum(in)/size(in,1));
handles.ClusterContrib=ClusterContrib;

sortedlist=cell(1,size(I,2));
for i=1:size(I,2);
    sortedlist(i)=strcat({'Cluster '},num2str(I(i)),{' - '},num2str(ClusterContrib(I(i),3)),{'%'});
end

if handles.ManualClusterCount==1;
    handles.cluster_sel.Value=[];
    handles.cluster_sel.String=sortedlist;
    handles.idx=double(in);
    handles.I=I;
    handles.Imod=I;
else
    p=handles.ManualClusterCount;
    idx=handles.idx+p*double(in);
    handles.idx=double(idx);
    handles.cluster_sel.Value=[];
    handles.cluster_sel.String=sortedlist;
    handles=HeatMap_CallbackManual(hObject, eventdata, handles);
    
end

guidata(hObject,handles);


% --- Executes on selection change in cluster_sel.
function cluster_sel_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cluster_sel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cluster_sel



% --- Executes during object creation, after setting all properties.
function cluster_sel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cluster_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearclusters.
function clearclusters_Callback(hObject, eventdata, handles)
% hObject    handle to clearclusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fieldremove={'ManualClusterCount','colorspec1','idx','ClusterMethod','ClusterContrib','num_clusters','HeatMapData','RowLabels','SizeCluster','I','Imod','Ifinal','line','thresholdcurrent','thresholdbook','threshold_count'};

for i=1:size(fieldremove,2);
    if isfield(handles,fieldremove{i})
        handles=rmfield(handles,fieldremove{i});
    end
end

handles.threshlist.String={};
handles.cluster_sel.String={};
handles.cluster_sel.Value=[];
handles.clusterplot.String={};
handles.clusterplot.Value=[];
handles.clusterfreq.String='';
handles.popupmenu1.Value=1;
handles.popupmenu2.Value=1;
Y=handles.Y;
scatter(handles.axes1,Y(:,1),Y(:,2),'filled');
handles.axes1.XTickLabel={};
handles.axes1.YTickLabel={};
guidata(hObject,handles);


% --- Executes on selection change in popupmenu1.
function handles=popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

handles=SortClusters(handles);
handles=ApplyCurrentThresh(handles);
handles=ClusterCut(handles);
PlotSelectClusters(handles);
guidata(hObject,handles);

function handles=SortClusters(handles)
    
    HeatMapData=handles.HeatMapData;
    ListC=[1:size(HeatMapData,1)];
    clear HeatMapData
    
    for j=ListC;
        clusterselect=handles.idx==j;
        SizeCluster(j)=sum(clusterselect);
        clusterselect2=handles.transy2(clusterselect,:);
        if size(clusterselect2,1)==1
            HeatMapData(j,:)=clusterselect2;
        else
            HeatMapData(j,:)=median(clusterselect2); 
        end
    end
    
    channelsel=handles.popupmenu1.Value;
    ClusterContrib=handles.ClusterContrib;
    channelselstring=handles.ChannelsAll{channelsel};
    valsort=strmatch(channelselstring,handles.ChannelsAll,'exact');
    I=handles.I;


    if handles.sortbutton.Value==0
        [B,I2]=sortrows(HeatMapData,-valsort);
    else
        [B,I2]=sortrows(HeatMapData,valsort);
    end

    I2=transpose(intersect(I2,I,'stable'));
    I=I2;

    sortedlist=cell(1,size(I,2));
    for i=1:size(I,2);
        sortedlist(i)=strcat({'Cluster '},num2str(I(i)),{' - '},num2str(ClusterContrib(I(i),3)),{'%'});
    end

    handles.cluster_sel.Value=[];    
    handles.cluster_sel.String=sortedlist;
    handles.Imod=I;

function handles=ApplyCurrentThresh(handles)

if isfield(handles,'thresholdbook')
    thresholdbook=handles.thresholdbook;
    HeatMapData=handles.HeatMapData;

    for i=1:size(thresholdbook,2);
        threshold_indx(i)=strmatch(thresholdbook(i).Channel,handles.ChannelsAll,'exact');
        threshold_dir{i}=thresholdbook(i).direction;
        threshold_val(i)=thresholdbook(i).threshold;
    end

    ListC=[1:size(HeatMapData,1)];
    clear HeatMapData
    thresh_cut_ind=ListC;
    
    for j=ListC;
            clusterselect=handles.idx==j;
            SizeCluster(j)=sum(clusterselect);
            clusterselect2=handles.transy2(clusterselect,:);
            if size(clusterselect2,1)==1
                HeatMapData(j,:)=clusterselect2;
            else
                HeatMapData(j,:)=median(clusterselect2); 
            end
    end
    
    for i=1:size(threshold_indx,2);
            thresh_cut=HeatMapData(:,threshold_indx(i));
            eval(['thresh_cut=thresh_cut' threshold_dir{i} 'threshold_val(i);'])
            thresh_cut_ind=intersect(ListC(thresh_cut),thresh_cut_ind);
    end

    I=handles.Imod;
    I=intersect(I,thresh_cut_ind,'stable');
    sortedlist=cell(1,size(I,2));
    for i=1:size(I,2);
        sortedlist(i)=strcat({'Cluster '},num2str(I(i)),{' - '},num2str(handles.ClusterContrib(I(i),3)),{'%'});
    end

    handles.cluster_sel.Value=[];
    handles.cluster_sel.String=sortedlist;
    handles.Imod=I;
end

function handles=ClusterCut(handles)
    SizeCluster=handles.SizeCluster;
    ClusterContrib=handles.ClusterContrib;
    I=handles.Imod;
    cut=str2num(handles.clusterfreq.String)/100;
    if isempty(cut);
        cut=0;
    end
    FreqCluster=SizeCluster./size(handles.y,1);
    Keep=(FreqCluster>cut).*[1:size(FreqCluster,2)];
    Keep(Keep==0)=[];
    I=intersect(I,Keep,'stable');
    handles.Imod=I;

    if ~isempty(I)
        for i=1:size(I,2);
        sortedlist(i)=strcat({'Cluster '},num2str(I(i)),{' - '},num2str(ClusterContrib(I(i),3)),{'%'});
        end
    end

    if ~exist('sortedlist')
        sortedlist={};
    end

    handles.cluster_sel.Value=[];
    handles.cluster_sel.String=sortedlist;
    PlotSelectClusters(handles);

    function PlotSelectClusters(handles)
    Imod=handles.Imod;
    I=handles.I;
    idx=handles.idx;
    colorspec1=handles.colorspec1;
    Y=handles.Y;
    
    if isfield(handles,'ManualClusterCount');
         Ybase=Y(idx==0,:);
         scatter(handles.axes1,Ybase(:,1),Ybase(:,2),'filled');
         hold(handles.axes1);
         Y=Y(idx>0,:);
         idx=idx(idx>0);
         
         selidx=ismember(idx,Imod);
         
         
         clear colorscheme
         for i=1:size(idx,1);
            colorscheme(i,:)=colorspec1(idx(i)).spec;
         end
         
       
        YNL=Y(~selidx,:);
        colorschemeNL=colorscheme(~selidx,:);
        scatter(handles.axes1,YNL(:,1),YNL(:,2),[],colorschemeNL,'filled','MarkerFaceAlpha',0.05);
        YHL=Y(selidx,:);
        colorschemeHL=colorscheme(selidx,:);
        scatter(handles.axes1,YHL(:,1),YHL(:,2),[],colorschemeHL,'filled');
        handles.axes1.XLim=handles.tsne_xlim;
        handles.axes1.YLim=handles.tsne_ylim;
        handles.axes1.XTickLabel={};
        handles.axes1.YTickLabel={};
        
        hold(handles.axes1);
        
    else

        clear colorscheme
        for i=1:size(idx,1);
            colorscheme(i,:)=colorspec1(idx(i)).spec;
        end

        selclusteridx=ismember(idx,Imod);
        YHL=Y(selclusteridx,:);
        colorschemeHL=colorscheme(selclusteridx,:);
        scatter(handles.axes1,YHL(:,1),YHL(:,2),[],colorschemeHL,'filled');
        hold(handles.axes1);
        YNL=Y(~selclusteridx,:);
        colorschemeNL=colorscheme(~selclusteridx,:);
        scatter(handles.axes1,YNL(:,1),YNL(:,2),[],colorschemeNL,'filled','MarkerFaceAlpha',0.05);
       
        handles.axes1.XTickLabel={};
        handles.axes1.YTickLabel={};
        handles.axes1.XLim=handles.tsne_xlim;
        handles.axes1.YLim=handles.tsne_ylim;
        hold(handles.axes1);
    
    end



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in HeatMapClusters.
function HeatMapClusters_Callback(hObject, eventdata, handles)
% hObject    handle to HeatMapClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'ManualClusterCount')
    ChannelsOut=handles.ChannelsOut;
    RowLabels=handles.RowLabels;
    HeatMapData=handles.HeatMapData;
    idx=handles.idx;
    y=handles.transy;
    
    
    select=handles.clusterplot.Value;
    if isfield(handles,'Ifinal')
        I=handles.Ifinal(select);
        if size(I,2)==2
            ChannelsOut=DetermineSigChannels(I,y,idx,ChannelsOut);
        end
        HeatMapData=HeatMapData(I,:);
        RowLabels=I;
    end
    
    hmobj=PlotHeatMap(HeatMapData,RowLabels,ChannelsOut);
    
    %PlotHeatMap(datasample(y,5000),[],ChannelsOut);
    
    
%     NewColumnLabels=hmobj.ColumnLabels;
%     NewRowLabels=hmobj.RowLabels;
%     
%     NewRowLabels=cellfun(@str2num,NewRowLabels);
%     %NewRowLabels=flip(NewRowLabels);
%     NewHeatMap=handles.HeatMapData(NewRowLabels,:);
%     
%     for i=1:size(NewColumnLabels,2);
%         temp=NewColumnLabels{i};
%         indx(i)=strmatch(temp,ChannelsOut,'exact');
%     end
%     
%     NewHeatMap=NewHeatMap(:,indx);
    
%     granularity=21;
%     colormap=redblue(granularity);
%     HMobj = HeatMap(NewHeatMap,'ColumnLabels',ChannelsOut,'RowLabels',RowLabels,'Standardize',1,'DisplayRange',3,'Colormap',colormap,'ColumnLabelsRotate',45);
%     
%     

% %     
%     indxfinal=[];
%     ClusterIndx=[];
%     n=1;
%     for i=transpose(NewRowLabels)
%         clusterselect=idx==i;
%         clusterselect2=y(clusterselect,:);
%         ClusterIndx=[ClusterIndx;n*ones(size(clusterselect2,1),1)];
%         n=n+1;
%         
%         if size(clusterselect2,1)==1
%             indxfinal=[indxfinal;find(clusterselect)];
%             continue
%         end
%         
%         indxorig=find(clusterselect);
% %         if size(clusterselect2,1)>50;
% %             [clusterselect2samp,idxsamp] = datasample(clusterselect2,50);
% %         else
% %             [clusterselect2samp,idxsamp] = datasample(clusterselect2,size(clusterselect2,1));
% %         end
% %         
% %         indxorig=indxorig(idxsamp); 
%         
%         Y=pdist(clusterselect2);
%         Z=linkage(Y);
%         leafOrder = optimalleaforder(Z,Y);
%         indxnew=indxorig(leafOrder);
%         indxfinal=[indxfinal;indxnew];
%     end
% % %     
%     ynewheat=y(indxfinal,:);
%     ynewheat=real(log10(ynewheat));
%     ynewheat=ynewheat(:,indx);
%     ynewheatz=zscore(ynewheat,0,1);
%         
%     
%     
%     n=n-1;
%     
%     figure;
%     for i=1:n
%         subplot(n,1,i);
%         granularity=21;
%         map=redblue(granularity);
%         colormap(map);
%         imagesc(ynewheatz(find(ClusterIndx==i),:));
%     end
%         
    
%     granularity=21;
%     colormap=redblue(granularity);
%     HMobj = HeatMap(ynewheatz,'ColumnLabels',NewColumnLabels,'Standardize',3,'DisplayRange',3,'Colormap',colormap,'ColumnLabelsRotate',45);
    
    
% %     
% %       granularity=21;
% %     map=redblue(granularity);
% %     figure
% %     colormap(map);
% %     imagesc(real(log10(ynewheat)));
% %     
%     
%     PlotExpandHeatMap(NewHeatMap,NewColumnLabels);
%     figure; heatmapcustom(ynewheat,NewColumnLabels,[],[], 'Colormap' , colormap, 'UseLogColormap',false,'ShowAllTicks',true,'ColorLevels',30)
    
    
    
else
    ChannelsOut=handles.ChannelsOut;
    idx=handles.idx;
    ClusterIter=[1:max(idx)];
    y=handles.transy;
    
    for j=ClusterIter
            clusterselect=idx==j;
            SizeCluster(j)=sum(clusterselect);
            clusterselect2=y(clusterselect,:);
            HeatMapData(j,:)=median(clusterselect2); 
            RowLabels{j}=strcat('Cluster ',num2str(j),' = ',num2str(100*(SizeCluster(j)/size(y,1))),'%');
    end
    

    PlotHeatMap(HeatMapData,RowLabels,ChannelsOut);
    handles.HeatMapData=HeatMapData;
    handles.RowLabels=RowLabels;
    ClusterContrib=tabulate(idx);
    if ClusterContrib(1,1)==0
        ClusterContrib=ClusterContrib(2:end,:);
    end
    handles.SizeCluster=SizeCluster;
    handles.ClusterContrib=ClusterContrib;
    handles.I=ClusterIter;
    handles.num_clusters=max(idx);
    guidata(hObject,handles);
    
    
end

% --- Executes on button press in HeatMapClusters.
function handles=HeatMap_CallbackManual(hObject, eventdata, handles)
% hObject    handle to HeatMapClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    ChannelsOut=handles.ChannelsOut;
    idx=handles.idx;
    ClusterIter=[1:max(idx)];
    y=handles.y;
    
    for j=ClusterIter
            clusterselect=idx==j;
            SizeCluster(j)=sum(clusterselect);
            clusterselect2=y(clusterselect,:);
            clusterselect2=log10(clusterselect2);
            HeatMapData(j,:)=real(median(clusterselect2)); 
            RowLabels{j}=strcat('Cluster ',num2str(j),' = ',num2str(100*(SizeCluster(j)/size(y,1))),'%');
    end
    
    handles.HeatMapData=HeatMapData;
    handles.RowLabels=RowLabels;
    ClusterContrib=tabulate(idx);
    if ClusterContrib(1,1)==0
        ClusterContrib=ClusterContrib(2:end,:);
    end
    handles.SizeCluster=SizeCluster;
    handles.ClusterContrib=ClusterContrib;
    handles.I=ClusterIter;
    handles.num_clusters=max(idx);
       


% --- Executes on button press in sortbutton.
function sortbutton_Callback(hObject, eventdata, handles)
% hObject    handle to sortbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sortbutton

if get(hObject,'Value')
    handles.sortbutton.String='Ascending';
else
    handles.sortbutton.String='Descending';
end

popupmenu1_Callback(hObject, eventdata, handles);


% --- Executes on button press in boxplotbutton.
function boxplotbutton_Callback(hObject, eventdata, handles)
% hObject    handle to boxplotbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%y=handles.y;
I=handles.Ifinal;
idx=handles.idx;
select=handles.clusterplot.Value;
I=I(select);
colorspec1=handles.colorspec1;
ChannelsOut=handles.ChannelsOut;
y=handles.transy;


if size(I,2)==2
    ChannelsOut=DetermineSigChannels(I,y,idx,ChannelsOut);
end


figure
selection=1;
clusterselect2=datasample(y,round(selection*size(y,1)));
%clusterselect2=real(log10(clusterselect2));
%clusterselect2=CleanLogData(clusterselect2);
for z=1:size(clusterselect2,2);
    r=normrnd(z,0.15,size(clusterselect2,1),1);
    dscatter(r,clusterselect2(:,z));
    hold on;
end

%n=2;
for j=I
    clusterselect=idx==j;
    clusterselect2=y(clusterselect,:);
    %clusterselect2=real(log10(clusterselect2));
    bp=boxplot(clusterselect2,'Colors',colorspec1(j).spec,'Symbol','','Whisker',0,'Notch','on'); 
    hold on
    for i = 1:size(bp,2), set(bp(:,i),'linewidth',3); end
    %n=n+1;
end

title('Phenotypic Characterization of Clusters')
xlabel('Marker');
ylabel('MFI (arcsinh)');
xtickangle(45);
xticks(1:size(y,2));
set(gca,'xticklabel',ChannelsOut);
set(gca,'TickLabelInterpreter','none');
hold off


% --- Executes on button press in hdflowplot.
function hdflowplot_Callback(hObject, eventdata, handles)
% hObject    handle to hdflowplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%y=handles.y;
y=handles.transy;
I=handles.Ifinal;
idx=handles.idx;
select=handles.clusterplot.Value;
I=I(select);
colorspec1=handles.colorspec1;
ChannelsOut=handles.ChannelsOut;
Num_Channels=size(y,2);

if size(I,2)==2
    ChannelsOut=DetermineSigChannels(I,y,idx,ChannelsOut);
end


    h=figure
    selection=1;
    %y=real(log10(y));
    %[y,keep]=CleanLogData(y);
    %idx=idx(keep);

    for z=1:size(y,2);
        r=normrnd(z,0.15,size(y,1),1);
        rkeep(:,z)=r;
        dscatter(r,y(:,z));
        hold on;

    end

    for j=I
        clusterselect=idx==j;
        clusterselect2=y(clusterselect,:);
        rkeep2=rkeep(clusterselect,:);
        for z=1:size(clusterselect2,2);
            scatter(rkeep2(:,z),clusterselect2(:,z),20,colorspec1(j).spec,'filled');
            hold on;
        end
    end


  title('Phenotypic Characterization of Clusters')
xlabel('Marker');
ylabel('MFI (arcsinh)');
xtickangle(45);
xticks(1:size(y,2));
set(gca,'xticklabel',ChannelsOut);
set(gca,'TickLabelInterpreter','none');
hold off

    function [data_clean,sel]=CleanLogData(data)
        [row,col]=find(data==-Inf);
        sel=setdiff([1:size(data,1)],row);
        data_clean=data(sel,:);
        
    function ChannelsOut=DetermineSigChannels(I,y,idx,ChannelsOut);
        
        cluster1=idx==I(1);
        cluster1=y(cluster1,:);
        
        cluster2=idx==I(2);
        cluster2=y(cluster2,:);
        
        n=1;
        for j=1:size(y,2)
            %[h,p]=ttest2(cluster1(:,j),cluster2(:,j));
            [p,h]=ranksum(cluster1(:,j),cluster2(:,j));
            pval(j)=p;   
        end
        
        
        Q = mafdr(pval,'BHFDR',true);
        
        for j=1:size(y,2)
             if Q(j)<0.05 && Q(j)>=0.01
                    temp=ChannelsOut{j};
                    temp=strcat(temp,'*');
                    ChannelsOut{j}=temp;
             elseif Q(j)<0.01 && Q(j)>=0.001
                 temp=ChannelsOut{j};
                 temp=strcat(temp,'**');
                 ChannelsOut{j}=temp;
             elseif Q(j)<0.001
                 temp=ChannelsOut{j};
                 temp=strcat(temp,'***');
                 ChannelsOut{j}=temp;
             end
        end

% --- Executes on button press in cvflowplot.
function cvflowplot_Callback(hObject, eventdata, handles)
% hObject    handle to cvflowplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

I=handles.Ifinal;
idx=handles.idx;
select=handles.clusterplot.Value;
I=I(select);
colorspec1=handles.colorspec1;
xscale=handles.xscale.Value;
yscale=handles.yscale.Value;

ChannelsAll=handles.ChannelsAll;
xplotnum=handles.xaxis.Value;
yplotnum=handles.yaxis.Value;
y=handles.y2;
ytrans=handles.transy2;


if xscale==1
     plotx=y(:,xplotnum);
elseif xscale==2
    plotx=log10forflow(y(:,xplotnum));
else
     plotx=ytrans(:,xplotnum);
%     objx=logicleTransform(max(y(:,xplotnum)),2,5,1);
%     plotx=objx.transform(y(:,xplotnum));
end

if yscale==1
    ploty=y(:,yplotnum);
elseif yscale==2
    ploty=log10forflow(y(:,yplotnum));
    
     %ploty=subplus(real(log10(y(:,yplotnum))));
else
    ploty=ytrans(:,yplotnum);
%     objy=logicleTransform(max(y(:,yplotnum)),2,5,1);
%     ploty=objy.transform(y(:,yplotnum));
end

figure
hold on
dscatter(plotx,ploty);

% if xscale==3
%     ax = gca;
%     ax.XTick = objx.Tick;
%     ax.XTickLabel = objx.TickLabel;
% end
% 
% if yscale==3
%     ax = gca;
%     ax.YTick = objy.Tick;
%     ax.YTickLabel = objy.TickLabel;
% end

  
    n=1;
    for j=I
        clusterselect=idx==j;
        xclusterplot=plotx(clusterselect);
        yclusterplot=ploty(clusterselect);
        scatter(xclusterplot,yclusterplot,10,colorspec1(j).spec,'filled','MarkerFaceAlpha',1);
    end
    


    title('Phenotypic Characterization of Clusters')
    xlabel(ChannelsAll(xplotnum));
    ylabel(ChannelsAll(yplotnum));


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

%h=uipanel('parent')
contents=cellstr(get(hObject,'String'));
Value=contents{get(hObject,'Value')};
gui_pop(handles,Value);

function gui_pop(handles,Value)
    positiongui=handles.figure1.Parent.PointerLocation;
    
    fh = figure('units','pixels',...
                  'units','normalized',...
                  'position',[0.25 0.25 .3 .6],...
                  'menubar','none',...
                  'name','Define Threshold',...
                  'numbertitle','off',...
                  'resize','off');
              
    guidata(fh,handles);
    
    handles.Value=Value;

    h=axes('Parent',fh,'Position',[0.05 0.1 0.55 0.8]);
    h.XTickLabel={};
    h.YTickLabel={};
    handles.h=h;

    sld = uicontrol('Parent',fh,'Style','slider',...
            'Min',0,'Max',1,'Value',0.5,...
            'units','normalized',...
            'Position', [0.65 0.1 0.05 0.8],...
            'Callback', @slider,...
            'Tag','slider1'); 
        
    handles.sld=sld;

    thresholdval=uicontrol('Parent',fh,'Style','edit',...
        'String','Threshold Value',...
        'units','normalized',...
        'Position',[.1,.03,0.2,0.05],...
        'Callback',@thresholdvaluetag,...
        'Tag','thresholdvalue');
    
    handles.thresholdval=thresholdval;

    yscaleloc=uicontrol('Style','pop',...
        'String',{'arcsinh';'linear'},...
        'units','normalized',...
        'Position',[.4,0.03,0.2,0.05],...
        'Tag','yscale',...
        'Callback',@yscaleloc);
    
    handles.yscaleloc=yscaleloc;
    
    addabove=uicontrol('Style','push',...
        'String','Add Above Threshold',...
        'units','normalized',...
        'Position',[0.75,.65,.2,.1],...
        'Tag','addabove',...
        'Callback',@addabove);
    
    handles.addabove=addabove;
    
    addbelow=uicontrol('Style','push',...
    'String','Add Below Threshold',...
    'units','normalized',...
    'Position',[.75,.45,0.2,0.1],...
    'Tag','addbelow',...
    'Callback',@addbelow);

    handles.addbelow=addbelow;

    title=uicontrol('Style','text',...
        'String',handles.popupmenu2.String{handles.popupmenu2.Value},...
        'units','normalized',...
        'Position',[.225,.925,0.2,0.05],...
        'FontSize',16);
    
    y=handles.y2;
    ytrans=handles.transy2;
    if handles.yscaleloc.Value==1
        %clusterselect2=subplus(log10(subplus(clusterselect2)));
        clusterselect2=ytrans;
    else
        clusterselect2=y;
    end
    r=normrnd(1,0.15,size(clusterselect2,1),1);
    index=strmatch(Value,handles.ChannelsAll,'exact');
    axes(h);
    dscatter(r,clusterselect2(:,index));
    h.XTickLabel={};
    h.XLim=[0,2];
    dataquery=clusterselect2(:,index);
    h.YLim=[prctile(dataquery,0),prctile(dataquery,100)];
    guidata(fh,handles);

    
function thresholdvaluetag(hObject, eventdata)
    % hObject    handle to thresholdvaluetag (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of thresholdvaluetag as text
    %        str2double(get(hObject,'String')) returns contents of thresholdvaluetag as a double
    handles=guidata(hObject);
    if isfield(handles,'line')
        delete(handles.line);
    end
    linepos=str2num(get(hObject,'String'));
    h=imline(handles.h,[-10 10],[linepos linepos]);
    handles.line=h;
    handles.sld.Value=linepos/5;
    guidata(hObject,handles);

function slider(hObject, eventdata)
    % hObject    handle to slider3 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'Value') returns position of slider
    %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    handles=guidata(hObject);
    if isfield(handles,'line')
        delete(handles.line);
    end
    index=strmatch(handles.Value,handles.ChannelsAll,'exact');
    
    if handles.yscaleloc.Value==1
        dataquery=handles.transy2(:,index);
        up=prctile(dataquery,100);
        down=prctile(dataquery,0);
        factor=(up-down);
        linepos=get(hObject,'Value')*factor+down;
    else
        dataquery=handles.y2(:,index);
        up=prctile(dataquery,100);
        down=prctile(dataquery,0);
        factor=(up-down);
        linepos=get(hObject,'Value')*factor+down;
    end
    
    h=imline(handles.h,[-10 10],[linepos linepos]);
    handles.line=h;
    hpos=getPosition(h);
    if handles.yscaleloc.Value==1
        handles.thresholdcurrent=hpos(1,2);
    else
        handles.thresholdcurrent=hpos(1,2);
    end
    handles.thresholdval.String=num2str(handles.thresholdcurrent);
    guidata(hObject,handles);

function yscaleloc(hObject,eventdata)
        handles=guidata(hObject);
        ytrans=handles.transy2;
        ylin=handles.y2;
        if get(hObject,'Value')==1
            clusterselect2=ytrans;
        else
            clusterselect2=ylin;
        end
        r=normrnd(1,0.15,size(clusterselect2,1),1);
        index=strmatch(handles.Value,handles.ChannelsAll,'exact');
        axes(handles.h);
        dscatter(r,clusterselect2(:,index));
        handles.h.XTickLabel={};
        handles.h.XLim=[0,2];
        dataquery=clusterselect2(:,index);
        handles.h.YLim=[prctile(dataquery,0),prctile(dataquery,100)];
        guidata(hObject,handles);
        
    function addabove(hObject,eventdata)
        handles=guidata(hObject);
        handles=SortClusters(handles);
        thresholdvalue=str2num(handles.thresholdval.String);
        channelselstring=handles.Value;

        if ~isfield(handles,'threshold_count')
            handles.thresholdbook(1).Channel=channelselstring;
            handles.thresholdbook(1).threshold=thresholdvalue;
            handles.thresholdbook(1).direction='>';
            handles.threshlist.String=strcat(channelselstring,{' , '},num2str(thresholdvalue),{' , '},'>');
            handles.threshold_count=1;
        else
            currentlist=handles.threshlist.String;
            count=handles.threshold_count;
            handles.thresholdbook(count+1).Channel=channelselstring;
            handles.thresholdbook(count+1).threshold=thresholdvalue;
            handles.thresholdbook(count+1).direction='>';
            currentlist=[currentlist;strcat(channelselstring,{' , '},num2str(thresholdvalue),{' , '},'>')];
            handles.threshlist.String=currentlist;
            handles.threshold_count=count+1;
        end
        handles=ApplyCurrentThresh(handles);
        handles=ClusterCut(handles);
        guidata(findobj('Tag','popupmenu2'),handles);
        closereq
        
    function addbelow(hObject,eventdata)
        handles=guidata(hObject);
        handles=SortClusters(handles);
        thresholdvalue=str2num(handles.thresholdval.String);
        channelselstring=handles.Value;

        if ~isfield(handles,'threshold_count')
            handles.thresholdbook(1).Channel=channelselstring;
            handles.thresholdbook(1).threshold=thresholdvalue;
            handles.thresholdbook(1).direction='<';
            handles.threshlist.String=strcat(channelselstring,{' , '},num2str(thresholdvalue),{' , '},'<');
            handles.threshold_count=1;
        else
            currentlist=handles.threshlist.String;
            count=handles.threshold_count;
            handles.thresholdbook(count+1).Channel=channelselstring;
            handles.thresholdbook(count+1).threshold=thresholdvalue;
            handles.thresholdbook(count+1).direction='<';
            currentlist=[currentlist;strcat(channelselstring,{' , '},num2str(thresholdvalue),{' , '},'<')];
            handles.threshlist.String=currentlist;
            handles.threshold_count=count+1;
        end
        handles=ApplyCurrentThresh(handles);
        handles=ClusterCut(handles);
        guidata(findobj('Tag','popupmenu2'),handles);
        closereq

       
    


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in threshlist.
function threshlist_Callback(hObject, eventdata, handles)
% hObject    handle to threshlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns threshlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from threshlist


% --- Executes during object creation, after setting all properties.
function threshlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearthresh.
function clearthresh_Callback(hObject, eventdata, handles)
% hObject    handle to clearthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'thresholdbook')
    handles=rmfield(handles,'thresholdbook');
    handles=rmfield(handles,'threshold_count');
end
handles.threshlist.Value=[1];
handles.threshlist.String=[];
handles.I=[1:handles.num_clusters];
guidata(hObject,handles);
popupmenu1_Callback(hObject, eventdata, handles)


function handles=clusterfreq_Callback(hObject, eventdata, handles)
% hObject    handle to clusterfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clusterfreq as text
%        str2double(get(hObject,'String')) returns contents of clusterfreq as a double

handles=SortClusters(handles);
handles=ApplyCurrentThresh(handles);
handles=ClusterCut(handles);
PlotSelectClusters(handles)
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function clusterfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in xaxis.
function xaxis_Callback(hObject, eventdata, handles)
% hObject    handle to xaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xaxis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xaxis


% --- Executes during object creation, after setting all properties.
function xaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yaxis.
function yaxis_Callback(hObject, eventdata, handles)
% hObject    handle to yaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns yaxis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yaxis


% --- Executes during object creation, after setting all properties.
function yaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in xscale.
function xscale_Callback(hObject, eventdata, handles)
% hObject    handle to xscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xscale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xscale


% --- Executes during object creation, after setting all properties.
function xscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yscale.
function yscale_Callback(hObject, eventdata, handles)
% hObject    handle to yscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns yscale contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yscale


% --- Executes during object creation, after setting all properties.
function yscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savetsneimg.
function savetsneimg_Callback(hObject, eventdata, handles)
% hObject    handle to savetsneimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path,filetype]=uiputfile({'*.bmp','BMP';'*.jpeg','JPEG';'*.png','PNG'},'Save Image As');
F=getframe(handles.axes1);
Image=frame2im(F);
imwrite(Image,strcat(path,file));


% --- Executes on selection change in clusterplot.
function clusterplot_Callback(hObject, eventdata, handles)
% hObject    handle to clusterplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clusterplot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clusterplot


% --- Executes during object creation, after setting all properties.
function clusterplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectclustersbutton.
function selectclustersbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectclustersbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sel=handles.cluster_sel.Value;
Imod=handles.Imod;
I=Imod(sel);
ClusterContrib=handles.ClusterContrib;

sortedlist=cell(1,size(I,2));
for i=1:size(sel,2);
    sortedlist(i)=strcat({'Cluster '},num2str(I(i)),{' - '},num2str(ClusterContrib(I(i),3)),{'%'});
end

currentlist=transpose(handles.clusterplot.String);
finallist=[currentlist,sortedlist];
finallist=unique(finallist,'stable');
handles.clusterplot.Value=[];
handles.clusterplot.String=finallist;
if isfield(handles,'Ifinal');
    handles.Ifinal=unique([handles.Ifinal,I],'stable');
else
    handles.Ifinal=I;
end
guidata(hObject,handles);


% --- Executes on button press in removecluster.
function removecluster_Callback(hObject, eventdata, handles)
% hObject    handle to removecluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sel=handles.clusterplot.Value;
currentlist=transpose(handles.clusterplot.String);
removelist=transpose(handles.clusterplot.String(sel));

newlist=setdiff(currentlist,removelist,'stable');
handles.clusterplot.Value=[];
handles.clusterplot.String=newlist;

handles.Ifinal=setdiff(handles.Ifinal,handles.Ifinal(sel),'stable');
guidata(hObject,handles);


% --- Executes on button press in cleandata.
function cleandata_Callback(hObject, eventdata, handles)
% hObject    handle to cleandata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.numall=handles.num;
    positiongui=handles.figure1.Parent.PointerLocation;
    
    fhgat = figure('units','pixels',...
                  'units','normalized',...
                  'position',[0.25 0.25 .3 .6],...
                  'menubar','none',...
                  'name','Gate Population',...
                  'numbertitle','off',...
                  'resize','off');
              
    guidata(fhgat,handles);
    num=handles.num;
    numall=handles.numall;
    ChannelsAll=handles.ChannelsAll;
    h=axes('Parent',fhgat,'Position',[0.25 0.15 0.7 0.7]);
    
    xoption=uicontrol('Parent',fhgat,'Style','pop',...
        'String',ChannelsAll,...
        'units','normalized',...
        'Position',[0.5,0.05,.2,.05],...
        'Value',1,...
        'Tag','xoption',...
        'Callback',@xoptionfunc1);
    
    yoption=uicontrol('Parent',fhgat,'Style','pop',...
        'String',ChannelsAll,...
        'units','normalized',...
        'Position',[0.02,0.5,0.2,0.05],...
        'Value',3,...
        'Tag','yoption',...
        'Callback',@yoptionfunc1);
    
    xscalegat=uicontrol('Parent',fhgat,'Style','pop',...
        'String',{'Linear','Log10','Arcsinh'},...
        'units','normalized',...
        'Position',[0.5,0.015,0.2,0.05],...
        'Value',1,...
        'Tag','xscalegat',...
        'Callback',@xscalegatfunc);
    
    yscalegat=uicontrol('Parent',fhgat,'Style','pop',...
        'String',{'Linear','Log10','Arcsinh'},...
        'units','normalized',...
        'Position',[0.02,0.465,0.2,0.05],...
        'Value',1,...
        'Tag','yscalegat',...
        'Callback',@yscalegatfunc);
    
    gatebutton=uicontrol('Parent',fhgat,'Style','push',...
        'String','Gate Population',...
        'units','normalized',...
        'Position',[0.1,.9,.3,.1],...
        'Tag','gatepop',...
        'Callback',@gatepop);
    
    donegate=uicontrol('Parent',fhgat,'Style','push',...
        'String','Apply Gates',...
        'units','normalized',...
        'Position',[0.4,0.9,0.3,0.1],...
        'Tag','donegate',...
        'Callback',@donegate);
        
    
    GatePlot(numall,num)
    guidata(fhgat,handles);
        
    
        function GatePlot(numall,num)
            xscale=findobj('Tag','xscalegat');
            yscale=findobj('Tag','yscalegat');
            xoption=findobj('Tag','xoption');
            yoption=findobj('Tag','yoption');

            if xscale.Value==1
                 plotx=num(:,xoption.Value);
            elseif xscale.Value==2
                plotx=log10forflow(num(:,xoption.Value));
            else
                plotx=asinh(num(:,xoption.Value));
%                 objx=logicleTransform(max(num(:,xoption.Value)),2,5,1);
%                 plotx=objx.transform(num(:,xoption.Value));
            end

            if yscale.Value==1
                ploty=num(:,yoption.Value);
            elseif yscale.Value==2
                 ploty=log10forflow(num(:,yoption.Value));
            else
                ploty=asinh(num(:,yoption.Value));
%                 objy=logicleTransform(max(num(:,yoption.Value)),2,5,1);
%                 ploty=objy.transform(num(:,yoption.Value));
            end
            
            dscatter(plotx,ploty);
            
%             if xscale.Value==3
%                 ax = gca;
%                 ax.XTick = objx.Tick;
%                 ax.XTickLabel = objx.TickLabel;
%             end
            
%             if yscale.Value==3
%                 ax = gca;
%                 ax.YTick = objy.Tick;
%                 ax.YTickLabel = objy.TickLabel;
%             end
                


        function xscalegatfunc(hObject,eventdata)
            handles=guidata(hObject);
            GatePlot(handles.numall,handles.num)

        function yscalegatfunc(hObject,eventdata)
            handles=guidata(hObject);
            GatePlot(handles.numall,handles.num)

        function gatepop(hObject,eventdata)
            handles=guidata(hObject);
            sel=selectdata('SelectionMode','Lasso','Verify','on');
            if ~isempty(sel);
                handles.num=handles.num(sel,:);
            end
            guidata(findobj('Tag','gatepop'),handles);
            GatePlot(handles.numall,handles.num);
            
        function xoptionfunc1(hObject,eventdata)
            handles=guidata(hObject);
            GatePlot(handles.numall,handles.num);
            
        function yoptionfunc1(hObject,eventdata)
            handles=guidata(hObject);
            GatePlot(handles.numall,handles.num);
            
        function donegate(hObject,eventdata)
            handles=guidata(hObject);
            guidata(findobj('Tag','cleandata'),handles);
            closereq
            
                
            


% --- Executes on button press in autocomp.
function autocomp_Callback(hObject, eventdata, handles)
% hObject    handle to autocomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(msgbox('Select Directory with Single Stains. Files should be labeled as Channel.fcs as indicated in the main name designation in the fcs file. Label unstained control as Unstained.fcs.'));
folder_name = uigetdir;
list=dir([strcat(folder_name)]);

        uiwait(msgbox('Select File to set Forward/SSC Gate'));

        [filename,folder]=uigetfile('*.fcs','Select file');
        filename=fullfile(folder,filename);

        [data,header]=fca_readfcs(filename);

        for j=1:size(header.par,2)
            temp=header.par(j).name;
            temp=strsplit(temp,'-A');
            ChannelsAll{j}=temp{1};
        end
        
        handles.data=data;
        handles.folder_name_single=folder_name;
        handles.list=list;
        fhgat = figure('units','normalized',... 
                       'position',[0.25 0.25 .3 .6],...
                      'menubar','none',...
                      'name','FSC/SSC Gating)',...
                      'numbertitle','off',...
                      'resize','off','Tag','fhgatfig');
                  
          guidata(fhgat,handles);

         h=axes('Parent',fhgat,'Position',[0.25 0.15 0.7 0.7]);

        xoption=uicontrol('Parent',fhgat,'Style','pop',...
            'String',ChannelsAll,...
            'units','normalized',...
            'Position',[0.5,0.05,.2,.05],...
            'Value',1,...
            'Tag','xoptioncomp',...
            'Callback',@xoptioncompfunc1);

        yoption=uicontrol('Parent',fhgat,'Style','pop',...
            'String',ChannelsAll,...
            'units','normalized',...
            'Position',[0.02,0.5,0.2,0.05],...
            'Value',3,...
            'Tag','yoptioncomp',...
            'Callback',@yoptioncompfunc1);

        xscalegat=uicontrol('Parent',fhgat,'Style','pop',...
            'String',{'Linear','Log10'},...
            'units','normalized',...
            'Position',[0.5,0.015,0.2,0.05],...
            'Value',1,...
            'Tag','xscalecompgat',...
            'Callback',@xscalegatcompfunc);

        yscalegat=uicontrol('Parent',fhgat,'Style','pop',...
            'String',{'Linear','Log10'},...
            'units','normalized',...
            'Position',[0.02,0.465,0.2,0.05],...
            'Value',1,...
            'Tag','yscalecompgat',...
            'Callback',@yscalegatcompfunc);

        gatebutton=uicontrol('Parent',fhgat,'Style','push',...
            'String','Gate Population',...
            'units','normalized',...
            'Position',[0.1,.9,.3,.1],...
            'Tag','gatepopcomp',...
            'Callback',@gatepopcomp);

        donegate=uicontrol('Parent',fhgat,'Style','push',...
            'String','Apply Gates',...
            'units','normalized',...
            'Position',[0.4,0.9,0.3,0.1],...
            'Tag','donegatecomp',...
            'Callback',@donegatecompfunc);

        GatePlotComp(data)
        
        


function GatePlotComp(data)
            xscale=findobj('Tag','xscalecompgat');
            yscale=findobj('Tag','yscalecompgat');
            xoption=findobj('Tag','xoptioncomp');
            yoption=findobj('Tag','yoptioncomp');

            if xscale.Value==1
                 func_callx='';
            else
                func_callx='log10';
            end

            if yscale.Value==1
                func_cally='';
            else
                func_cally='log10'; 
            end
            
            
            
            eval(['dscatter(' func_callx '(data(:,xoption.Value)),' func_cally '(data(:,yoption.Value)));'])


        function xscalegatcompfunc(hObject,eventdata)
            handles=guidata(hObject);
            GatePlotComp(handles.data)


        function yscalegatcompfunc(hObject,eventdata)
            handles=guidata(hObject);
            GatePlotComp(handles.data);

        function gatepopcomp(hObject,eventdata)
            handles=guidata(hObject);
            [sel,xsel,ysel]=selectdata('SelectionMode','Lasso','Verify','on');
            handles.xsel=xsel;
            handles.ysel=ysel;
            guidata(findobj('Tag','gatepopcomp'),handles);
            GatePlotComp(handles.data);
            donegatecompfunc(hObject,eventdata)

            
        function xoptioncompfunc1(hObject,eventdata)
            handles=guidata(hObject);
            GatePlotComp(handles.data);
 
        function yoptioncompfunc1(hObject,eventdata)
            handles=guidata(hObject);
            GatePlotComp(handles.data);

            
        function donegatecompfunc(hObject,eventdata)
            handles=guidata(hObject);
             xscale=findobj('Tag','xscalecompgat');
            yscale=findobj('Tag','yscalecompgat');
            xoption=findobj('Tag','xoptioncomp');
            yoption=findobj('Tag','yoptioncomp');
            
            if xscale.Value==2
                xsel=10.^handles.xsel;
            else
                xsel=handles.xsel;
            end
            if yscale.Value==2
                ysel=10.^handles.ysel;
            else
                ysel=handles.ysel;
            end
            
            xoption_val=xoption.Value;
            yoption_val=yoption.Value;
            xscale_val=xscale.Value;
            yscale_val=yscale.Value
            
            
            closereq
            if isfield(handles,'data')
                handles=rmfield(handles,'data')
            end
            
            
            k=boundary(xsel,ysel);
            boundx=xsel(k);
            boundy=ysel(k);

            [data_unstained,header_unstained]=fca_readfcs(fullfile(handles.folder_name_single,'Unstained.fcs'));
            in=inpolygon(data_unstained(:,xoption_val),data_unstained(:,yoption_val),boundx,boundy);
            data_unstained=data_unstained(in,:);
            
            list=handles.list;

            for i=1:size(list,1)
                filename=list(i).name;
                [pathstr,name,ext] = fileparts(filename);
               
                if strcmp(ext,'.fcs')
                    if ~strcmp(filename,'Unstained.fcs')
                         filename2=fullfile(handles.folder_name_single,filename);
                        [data,header]=fca_readfcs(filename2);
                        in=inpolygon(data(:,xoption_val),data(:,yoption_val),boundx,boundy);
                        data=data(in,:);

                        for j=1:size(header.par,2)
                            temp=header.par(j).name;
                            temp=strsplit(temp,'-A');
                            ChannelsAll{j}=temp{1};
                        end

                        filenamequery=strsplit(filename,'.');
                        filenamequery=filenamequery{1};
                        loc=strmatch(filenamequery,ChannelsAll,'exact');
                        data_out=data(:,loc);
                        data_out_un=data_unstained(:,loc);
                        cut=prctile(data_out_un,99);

  
                        num_clusters=50;
                        idx=kmeans(data_out,num_clusters);

                        results=tabulate(idx);
                        for j=1:num_clusters
                            clusterselect=idx==j;
                            datacluster=data(clusterselect,:);
                            int_cluster=min(datacluster(:,loc));
                            results(j,4)=int_cluster;
                        end

                        keep=results(:,4)>cut;
                        results=results(keep,:);

                        clusterselect=ismember(idx,results(:,1));
                        datacluster=data(clusterselect,:);
                         Data(loc).pop=datacluster;
                    end
                end
            end
            
            n=1;
            for i=1:size(Data,2);
                if ~isempty(Data(i).pop)
                    comp_ind(n)=i;
                    n=n+1;
                end
            end

            num_param=size(comp_ind,2);

            a=mean(data_unstained(:,comp_ind));
            
            clear k
            for i=1:num_param
                for j=1:num_param
                    top=subplus(sum(Data(comp_ind(i)).pop(:,comp_ind(j)))-size(Data(comp_ind(i)).pop,1)*a(j));
                    bottom=subplus(sum(sum(Data(comp_ind(i)).pop(:,comp_ind))-size(Data(comp_ind(i)).pop,1).*a));
                    k(j,i)=top/bottom;
                end
            end
        
        DataSample=handles.num2;    
        O=transpose(DataSample(:,comp_ind));
        A=transpose(a);
        OS=O-A;
        S=inv(k)*(O-A);
        DataSample2=DataSample;
        DataSample2(:,comp_ind)=transpose(S);
        
%         
%         i=12;
%         n=1;
%         comp=[7:9];
%         figure
%         for j=comp
%             subplot(size(comp,2),2,n)
%             dscatter(transpose(log10(OS(i-6,:))),transpose(log10(OS(j-6,:))))
%             xlabel(ChannelsAll{i});
%             ylabel(ChannelsAll{j})
%             n=n+1;
%             subplot(size(comp,2),2,n)
%             dscatter(log10(DataSample2(:,i)),log10(DataSample2(:,j)))
%             xlabel(header.par(i).name2);
%             n=n+1;
%         end
        
        handles.num=DataSample2;
        handles=rmfield(handles,'folder_name_single');
        handles=rmfield(handles,'list');
        handles=rmfield(handles,'xsel');
        handles=rmfield(handles,'ysel');
        guidata(findobj('Tag','autocomp'),handles);
        msgbox('Compensation Completed');
        
            
         

% --- Executes on button press in distinguish_clusters.
function distinguish_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to distinguish_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


I=handles.Ifinal;
select=handles.clusterplot.Value;
I=I(select);

HeatMapData=handles.HeatMapData;
HeatMapData=HeatMapData(I,:);
Z=zscore(HeatMapData);
Zhi=Z>1;
Zlo=Z<1;


% --- Executes on selection change in heatmaptsne.
function heatmaptsne_Callback(hObject, eventdata, handles)
% hObject    handle to heatmaptsne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns heatmaptsne contents as cell array
%        contents{get(hObject,'Value')} returns selected item from heatmaptsne

contents=cellstr(get(hObject,'String'));
sel=get(hObject,'Value');
Y=handles.Y;
y2=handles.transy2;
channel=y2(:,sel);

cutofftop=prctile(channel,100);
cutoffbottom=prctile(channel,0);
replace=(channel>cutofftop);
channel(replace)=cutofftop;
replace=channel<cutoffbottom;
channel(replace)=cutoffbottom;


figure('Name',contents{sel},'NumberTitle','off');
scatter(Y(:,1),Y(:,2),5,channel,'filled');
colormap('jet');
title(contents(sel),'FontSize',20,'Interpreter','none')


% --- Executes during object creation, after setting all properties.
function heatmaptsne_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heatmaptsne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveworkspace.
function saveworkspace_Callback(hObject, eventdata, handles)
% hObject    handle to saveworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[file,path,filetype]=uiputfile({'*.mat','MAT'});
if ~path==0
    hbox=msgbox('Saving Workspace...');
    save(strcat(path,file),'handles');
    close(hbox);
else
    return
end


% --- Executes on button press in heatmapind.
function heatmapind_Callback(hObject, eventdata, handles)
% hObject    handle to heatmapind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ChannelsOut=handles.ChannelsAll(handles.channel_select);
y=handles.transy;
%HeatMapData=handles.HeatMapData;
%idx=handles.idx;
%RowLabels=[];
%hmobj=PlotHeatMap(HeatMapData,RowLabels,ChannelsOut);
    
    
    
%     NewColumnLabels=hmobj.ColumnLabels;
%     NewRowLabels=hmobj.RowLabels;
% %     
%     NewRowLabels=cellfun(@str2num,NewRowLabels);
%     NewHeatMap=HeatMapData(NewRowLabels,:);
% %     
%     for i=1:size(NewColumnLabels,2);
%         temp=NewColumnLabels{i};
%         indx(i)=strmatch(temp,ChannelsOut,'exact');
%     end
% %     
%     NewHeatMap=NewHeatMap(:,indx);
%     %HeatMap(NewHeatMap,'ColumnLabels',NewColumnLabels,'Standardize',1,'DisplayRange',3,'Colormap',colormap)
% % %     
%     indxfinal=[];
%     for i=transpose(NewRowLabels)
%         clusterselect=idx==i;
%         clusterselect2=y(clusterselect,:);
%         
%         if size(clusterselect2,1)==1
%             indxfinal=[indxfinal;find(clusterselect)];
%             continue
%         end
%         
%         indxorig=find(clusterselect);
%         if size(clusterselect2,1)>50;
%             [clusterselect2samp,idxsamp] = datasample(clusterselect2,50);
%         else
%             [clusterselect2samp,idxsamp] = datasample(clusterselect2,size(clusterselect2,1));
%         end
%         
%         indxorig=indxorig(idxsamp); 
%         
%         Y=pdist(clusterselect2samp);
%         Z=linkage(Y);
%         leafOrder = optimalleaforder(Z,Y);
%         indxnew=indxorig(leafOrder);
%         indxfinal=[indxfinal;indxnew];
%     end
% % %     
%      ynewheat=y(indxfinal,:);
%      ynewheat=real(log10(ynewheat));
% % %     
%         granularity=21;
%         colormap=redblue(granularity);
%         HeatMap(ynewheat,'ColumnLabels',NewColumnLabels,'Standardize',1,'DisplayRange',3,'Colormap',colormap)
%     

hbox=msgbox('Calculating HeatMap...');
PlotHeatMap(y,[],ChannelsOut);
close(hbox);


% --- Executes on button press in gatetsne.
function gatetsne_Callback(hObject, eventdata, handles)
% hObject    handle to gatetsne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dcm_obj=datacursormode(handles.axes1.Parent);
set(dcm_obj,'Enable','off');

sel=selectdata('Axes',handles.axes1,'SelectionMode','Lasso','Verify','on');
if iscell(sel)
    sel=sel{end};
end

if isfield(handles,'Imod');
    Imod=handles.Imod;
    idx=handles.idx;
    selclusteridx=ismember(idx,Imod);
    selclusteridx=find(selclusteridx);
    sel2=selclusteridx(sel);
    sel=sel2;
end

handles.y2=handles.y2(sel,:);
handles.y=handles.y(sel,:);
handles.transy2=handles.transy2(sel,:);
handles.transy=handles.transy(sel,:);

handles.num_samples=size(handles.y,1);
handles.numsamples.String=num2str(handles.num_samples);
handles.perc_file.String=num2str(100*(handles.num_samples/size(handles.num,1)));


normalizetsne=1;

hbox=msgbox('Running t-SNE Analysis');
Y=tsne(handles.transy,'Standardize',normalizetsne);
close(hbox);
handles.Y=Y;

scatter(handles.axes1,Y(:,1),Y(:,2),'filled');
handles.axes1.XTickLabel={};
handles.axes1.YTickLabel={};
handles.axes1.Title.String='TSNE';
handles.tsne_xlim=handles.axes1.XLim;
handles.tsne_ylim=handles.axes1.YLim;

guidata(hObject,handles);


% --- Executes on button press in loadworkspace.
function loadworkspace_Callback(hObject, eventdata, handles)
% hObject    handle to loadworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path,filetype]=uigetfile({'*.mat','MAT'});
if ~path==0
load(strcat(path,file),'handles');
currentobjects=findall(0);
else
    return
end

n=1;
for i=1:size(currentobjects,1);
    try
    nametemp=currentobjects(i).Name;
    indx(n)=i;
    n=n+1;
    catch
        continue
    end
end

close(currentobjects(indx(end)));






function eventperfile_Callback(hObject, eventdata, handles)
% hObject    handle to eventperfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eventperfile as text
%        str2double(get(hObject,'String')) returns contents of eventperfile as a double


% --- Executes during object creation, after setting all properties.
function eventperfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventperfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in inst_type.
function inst_type_Callback(hObject, eventdata, handles)
% hObject    handle to inst_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns inst_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from inst_type


% --- Executes during object creation, after setting all properties.
function inst_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inst_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
