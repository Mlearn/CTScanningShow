function ScanGeoPlotFigNew()
% SCANGEOPLOTFIGNEW CT??????????????
% Created by An Mou
% Oct. 13, 2016

global handles

GeoParaCreate(); % ????????
UI_Build();      % ????UI????
GeoParaConectUI(); %??????????????
UIConectGeoPara(); %??????????????
UIValueInit();     %??????UI??????????
TotalPlot();       %????????

handles.AutoPlayhandle = timer('Name','CircleTimer', 'TimerFcn', @AutoPlayTask, 'Period',0.2, 'ExecutionMode','fixedspacing');
handles.AutoPlayStatus = 'Original';
end
function UI_Build()
global handles
scrsz = get(0,'ScreenSize');

main_fig_size_position = [scrsz(1)+floor(0.3*scrsz(3)), scrsz(2)+floor(0.3*scrsz(4)), floor(0.5*scrsz(3)), floor(0.55*scrsz(4))];
handles.main_fig = figure('Tag','main_fig','Position',main_fig_size_position, 'Resize','on');
handles.ShowPanel = axes('visible','on','Position',[0 0 0.7 1],'XTick',[],'YTick',[],'ZTick',[],'Children',[],'Box','on');%,'CameraPosition',[1,0.1,0.05]);
handles.ActPanel = axes('Visible','off','Position',[0.7 0 0.19 1]);
handles.UI.AngViewSD.handle=uicontrol(handles.main_fig,'style','slider','Min',0,'Max',10,'Value',0,'SliderStep',[0.1,0.1],'Units','normalized',...
    'Position',[0.97 0 0.03 1],'callback',{@Scan_Callback,'AngViewSD'});

xL1 = 0.71; xW1 = 0.12; xL2 = xL1+xW1; xW2 = 0.08; h = 0.04; yL = 0.9:-h-0.02:0.2;
PoName = [repmat(xL1,[length(yL),1]), yL',repmat([xW1,h],[length(yL),1])];
PoValu = [repmat(xL2,[length(yL),1]), yL',repmat([xW2,h],[length(yL),1])];
%====Geometry Parameter
handles.UI.GeometryName.handle=uicontrol(handles.main_fig,'Style','text','String','Geometry',...
    'HorizontalAlignment','center','FontWeight','Bold','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(1,:)+[0 0 xW2 0]);
% ==??????????
handles.UI.SliceNumName.handle=uicontrol(handles.main_fig,'Style','text','String','SliceNum:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(2,:));
handles.UI.SliceNumEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','32',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(2,:),'callback',{@Edit_GeometryReset_Callback, 'SliceNumEdit'});
% ==??????????????
handles.UI.SliceWidthName.handle=uicontrol(handles.main_fig,'Style','text','String','SliceWidth:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(3,:));
handles.UI.SliceWidthEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','0.625',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(3,:),'callback',{@Edit_GeometryReset_Callback, 'SliceWidthEdit'});
handles.UI.SliceWidthUnit.handle=uicontrol(handles.main_fig,'Style','text','String','mm',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(3,:)+[xW2 0 -xW2+0.06 0]);
% ==????
handles.UI.PitchName.handle=uicontrol(handles.main_fig,'Style','text','String','Pitch:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(4,:));
handles.UI.PitchEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','0.8',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(4,:),'callback',{@Edit_GeometryReset_Callback, 'PitchEdit'});

%====Scan Parameter
handles.UI.ScanName.handle=uicontrol(handles.main_fig,'Style','text','String','Scan',...
    'HorizontalAlignment','center','FontWeight','Bold','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(5,:)+[0 0 xW2 0]);
% ==????????????????
handles.UI.AngViewMinName.handle=uicontrol(handles.main_fig,'Style','text','String','ViewMin:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(6,:));
handles.UI.AngViewMinEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','0',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(6,:),'callback',{@Edit_GeometryReset_Callback, 'AngViewMinEdit'});
handles.UI.AngViewMinUnit.handle=uicontrol(handles.main_fig,'Style','text','String','??',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(6,:)+[xW2 0 -xW2+0.04 0]);
% ==????????????????
handles.UI.AngViewMaxName.handle=uicontrol(handles.main_fig,'Style','text','String','ViewMax:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(7,:));
handles.UI.AngViewMaxEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','5',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(7,:),'callback',{@Edit_GeometryReset_Callback, 'AngViewMaxEdit'});
handles.UI.AngViewMaxUnit.handle=uicontrol(handles.main_fig,'Style','text','String','??',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(7,:)+[xW2 0 -xW2+0.04 0]);
% ==????????????
handles.UI.AngViewStepName.handle=uicontrol(handles.main_fig,'Style','text','String','ViewStep:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(8,:));
handles.UI.AngViewStepEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','0.01',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(8,:),'callback',{@Edit_GeometryReset_Callback, 'AngViewStepEdit'});
handles.UI.AngViewStepUnit.handle=uicontrol(handles.main_fig,'Style','text','String','??',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(8,:)+[xW2 0 -xW2+0.04 0]);
% ==????????
handles.UI.AngViewName.handle=uicontrol(handles.main_fig,'Style','text','String','ProjView:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(9,:));
handles.UI.AngViewEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','0',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(9,:),'callback',{@Scan_Callback, 'AngViewEdit'});
handles.UI.AngViewUnit.handle=uicontrol(handles.main_fig,'Style','text','String','??',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(9,:)+[xW2 0 -xW2+0.04 0]);

% ==
handles.UI.AutoPlayButton.Handle = uicontrol('Style', 'pushbutton', 'String', 'Play', 'Units','normalized',...
        'Position', [0.72 0.04 0.08 0.06], 'Callback', {@Button_AutoPlay_Callback, 'Play'});
handles.UI.AutoPlayButton.Handle = uicontrol('Style', 'pushbutton', 'String', 'Pause', 'Units','normalized',...
        'Position', [0.8 0.04 0.08 0.06], 'Callback', {@Button_AutoPlay_Callback, 'Pause'});
handles.UI.AutoPlayButton.Handle = uicontrol('Style', 'pushbutton', 'String', 'Reset', 'Units','normalized',...
        'Position', [0.88 0.04 0.08 0.06], 'Callback', {@Button_AutoPlay_Callback, 'Reset'});
end

function GeoParaCreate()
global handles
handles.GeoPara.sourceToIso.Value = 570;
handles.GeoPara.Pitch.Value = 0.8;
handles.GeoPara.SliceWidth.Value = 0.625;
handles.GeoPara.SliceNum.Value = 32;
handles.GeoPara.isoToDetector.Value = 445;
handles.GeoPara.Zanchor.Value = 0;
handles.GeoPara.AngViewMin.Value = 0;
handles.GeoPara.AngViewMax.Value = 5;
handles.GeoPara.AngViewStep.Value = 0.1;
handles.GeoPara.AngViewValue.Value = 0;
handles.GeoPara.SFOV.Value = 500; %mm

Samp = 16;
handles.GeoPara.ChannelNumPar.Value = 1665/Samp;
handles.GeoPara.MidChannelPar.Value = 8.265021361983943e+02/Samp;
handles.GeoPara.ChannelParSpace.Value = 0.305497508156251*Samp;
end
function GeoParaConectUI()
global handles
handles.GeoPara.sourceToIso.relateUI = {};
handles.GeoPara.Pitch.relateUI = {'PitchEdit'};
handles.GeoPara.SliceWidth.relateUI = {'SliceWidthEdit'};
handles.GeoPara.SliceNum.relateUI = {'SliceNumEdit'};
handles.GeoPara.isoToDetector.relateUI = {};
handles.GeoPara.Zanchor.relateUI = {};
handles.GeoPara.AngViewMin.relateUI = {'AngViewMinEdit','AngViewSD'};
handles.GeoPara.AngViewMax.relateUI = {'AngViewMaxEdit','AngViewSD'};
handles.GeoPara.AngViewStep.relateUI = {'AngViewStepEdit'};
handles.GeoPara.AngViewValue.relateUI = {'AngViewEdit','AngViewSD'};
handles.GeoPara.SFOV.relateUI = {};

handles.GeoPara.ChannelNumPar.relateUI = {};
handles.GeoPara.MidChannelPar.relateUI = {};
handles.GeoPara.ChannelParSpace.relateUI = {};
end
function UIConectGeoPara()
global handles
handles.UI.PitchEdit.relatePara = {'Pitch'};
handles.UI.SliceWidthEdit.relatePara = {'SliceWidth'};
handles.UI.SliceNumEdit.relatePara = {'SliceNum'};
handles.UI.AngViewMinEdit.relatePara = {'AngViewMin'};
handles.UI.AngViewMaxEdit.relatePara = {'AngViewMax'};
handles.UI.AngViewStepEdit.relatePara = {'AngViewStep'};
handles.UI.AngViewEdit.relatePara = {'AngViewValue'};
handles.UI.AngViewSD.relatePara = {'AngViewValue'};
end
function GeoParaValueSet(ParaName, Value)
global handles
if strcmp('AngViewMin', ParaName)
    handles.GeoPara.AngViewMin.Value = min(Value, handles.GeoPara.AngViewMax.Value);
    handles.GeoPara.AngViewValue.Value = max(handles.GeoPara.AngViewValue.Value, handles.GeoPara.AngViewMin.Value);
elseif strcmp('AngViewMax', ParaName)
    handles.GeoPara.AngViewMax.Value = max(Value, handles.GeoPara.AngViewMin.Value);
    handles.GeoPara.AngViewValue.Value = min(handles.GeoPara.AngViewValue.Value, handles.GeoPara.AngViewMax.Value);
elseif strcmp('AngViewValue', ParaName)
    handles.GeoPara.AngViewValue.Value = min(Value, handles.GeoPara.AngViewMax.Value);
    handles.GeoPara.AngViewValue.Value = max(Value, handles.GeoPara.AngViewMin.Value);
else
    handles.GeoPara.(ParaName).Value = Value;
end
end
function SDUISet(UIName)
global handles
MinValue = handles.GeoPara.AngViewMin.Value;
MaxValue = handles.GeoPara.AngViewMax.Value;
Value = handles.GeoPara.AngViewValue.Value;
set(handles.UI.(UIName).handle,'Min',MinValue, 'Max', MaxValue, 'Value', Value);
end
function ParaName = GeoParaUpdateFromUI(UIName)
% ????????????UI??????????????????????????
global handles
ParaName = [];
if isfield(handles.UI.(UIName),'relatePara')
    ParaName = handles.UI.(UIName).relatePara{1};
    if strcmp('edit', handles.UI.(UIName).handle.Style) % ????UI??????????????
        GeoParaValueSet(ParaName, str2double(get(handles.UI.(UIName).handle,'string')));
    end
    if strcmp('slider', handles.UI.(UIName).handle.Style) % ????UI??????
        GeoParaValueSet(ParaName, get(handles.UI.(UIName).handle,'Value'));
    end
end
end
function UIUpdateFromGeoPara(ParaName)
% ??????GeoPara??????????????????????
global handles
if ~isempty(ParaName)
    if isfield(handles.GeoPara.(ParaName), 'relateUI')
        UINames = handles.GeoPara.(ParaName).relateUI;
        for i = 1:length(UINames)
            UIName = UINames{i};
            UIStyle = handles.UI.(UIName).handle.Style;
            if strcmp('slider', UIStyle)
                SDUISet(UIName);
            end
            if strcmp('edit', UIStyle)
                set(handles.UI.(UIName).handle, 'String', num2str(handles.GeoPara.(ParaName).Value));
            end
        end
    end
end
end
function UIValueInit()
global handles
ParaNames = fieldnames(handles.GeoPara);
for i = 1:length(ParaNames)
    ParaName = ParaNames{i};
    UIUpdateFromGeoPara(ParaName);
end
end

function Button_AutoPlay_Callback(hObject, callbackdata, Sig)
global handles
if strcmp('Play', Sig) && ~strcmp('Play', handles.AutoPlayStatus)
    start(handles.AutoPlayhandle);
    handles.AutoPlayStatus = 'Play';
elseif strcmp('Pause', Sig) && ~strcmp('Pause', handles.AutoPlayStatus)
    stop(handles.AutoPlayhandle);
    handles.AutoPlayStatus = 'Pause';
elseif strcmp('Reset', Sig)
    if strcmp('Play', handles.AutoPlayStatus)
    stop(handles.AutoPlayhandle);
    end
    GeoParaCreate(); % ????????
    UIValueInit();     %??????UI??????????
    TotalPlot();       %????????
    handles.AutoPlayStatus = 'Original';
end
end
function Scan_Callback(hObject, callbackdata, UIName)
global handles
ParaName = GeoParaUpdateFromUI(UIName);
UIUpdateFromGeoPara(ParaName);
DetecPlotRefresh();
end
function Edit_GeometryReset_Callback(hObject, callbackdata, UIName)
%The function definition must define two input arguments, hObject and callbackdata, 
%followed by any additional arguments the function uses. Handle Graphics automatically 
%passes the hObject and callbackdata arguments when it calls the function. If you define 
%additional input arguments, then the values you pass to the callback must exist in the 
%workspace when the end user triggers the callback.
global handles
ParaName = GeoParaUpdateFromUI(UIName);
UIUpdateFromGeoPara(ParaName);
TotalPlot();
end

function TotalPlot()
global handles
GeoParameter = handles.GeoPara;
cla(handles.ShowPanel);
TrajectoryPlot(handles.ShowPanel, GeoParameter);

fCouchDistPerView = GeoParameter.SliceWidth.Value*GeoParameter.SliceNum.Value*GeoParameter.Pitch.Value/2/pi;
AngleArray = GeoParameter.AngViewMin.Value:GeoParameter.AngViewStep.Value:GeoParameter.AngViewMax.Value;
AngleCenter = pi*(AngleArray(1)+AngleArray(end))/2;
PlanePlot(handles.ShowPanel, fCouchDistPerView*AngleCenter, GeoParameter);
handles.Detect = DetecPlot(handles.ShowPanel, GeoParameter);
drawnow;
end
function h = TrajectoryPlot(ParentH, GeoParameter)
% Geometry parameter
sourceToIso = GeoParameter.sourceToIso.Value;
Pitch = GeoParameter.Pitch.Value;
SliceWidth = GeoParameter.SliceWidth.Value;
SliceNum = GeoParameter.SliceNum.Value;
DetectorZ = SliceWidth*SliceNum;
Zanchor = GeoParameter.Zanchor.Value;
fCouchDistPerView = DetectorZ*Pitch/2/pi;

AngleArray = pi*(GeoParameter.AngViewMin.Value:0.01:GeoParameter.AngViewMax.Value);

Xarray = sourceToIso.*sin(AngleArray);
Yarray = sourceToIso.*cos(AngleArray);
Zarray = fCouchDistPerView*AngleArray+Zanchor;
axes(ParentH); hold on
h = plot3(Xarray,Yarray,Zarray,'r-');

end
function h = PlanePlot(ParentH, Zposi, GeoParameter)

axes(ParentH); hold on
%xy = -GeoParameter.SFOV.Value/2:GeoParameter.SFOV.Value/2;
ImgZ = Zposi;
%surf(xy,xy, ones(length(xy),length(xy))*ImgZ); % ????????????,????Parent????????????????????
hb = GeoParameter.SFOV.Value/2;
% vv = [-hb,-hb,ImgZ; hb,-hb,ImgZ;hb,hb,ImgZ;-hb,hb,ImgZ];
% patch('vertices', vv, 'faces', [1,2,3,4], 'FaceColor', 'y','FaceAlpha',0.5);
% x = 250*sin(0:0.03:2*pi);
% y = 250*cos(0:0.03:2*pi);
% z = ones(length(x))*ImgZ;
% h = plot3(x,y,z, '*g');
h = rectangle('Position',[-hb,-hb,2*hb,2*hb],'Curvature',[1,1],  'FaceColor','r');
end
function h = DetecPlot(parentH, GeoParameter)
Angle = pi*GeoParameter.AngViewValue.Value;
Pitch = GeoParameter.Pitch.Value;
sourceToIso = GeoParameter.sourceToIso.Value;
isoToDetector = GeoParameter.isoToDetector.Value;
SliceWidth = GeoParameter.SliceWidth.Value;
SliceNum = GeoParameter.SliceNum.Value;
DetectorZ = SliceWidth*SliceNum;
DetectorZArray = (sourceToIso+isoToDetector)/sourceToIso*((-(SliceNum-1)/2:(SliceNum-1)/2)*SliceWidth)';
Zanchor = GeoParameter.Zanchor.Value;
ChannelNumPar = GeoParameter.ChannelNumPar.Value;
MidChannelPar = GeoParameter.MidChannelPar.Value;
ChannelParSpace = GeoParameter.ChannelParSpace.Value;

fCouchDistPerView = DetectorZ*Pitch/2/pi;

ChannelParArrayX = ChannelParSpace*((1:ChannelNumPar)-MidChannelPar);
SourceArrayX = ChannelParArrayX;
SourceAngleArray = asin(SourceArrayX./sourceToIso);
SourceArrayY = sqrt(sourceToIso.^2-SourceArrayX.^2);
ChannelParArrayY = SourceArrayY-sourceToIso-isoToDetector;
ZShiftArray = SourceAngleArray .* fCouchDistPerView;
SourceArrayZ = ZShiftArray;
% plot3(SourceArrayX, SourceArrayY, SourceArrayZ);

DetectorArrayX = SourceArrayX;
DetectorArrayY = ChannelParArrayY;
DetectorArrayZ = ZShiftArray;

SourceArrayXr =  SourceArrayX.*cos(Angle) + SourceArrayY.*sin(Angle);
SourceArrayYr = -SourceArrayX.*sin(Angle) + SourceArrayY.*cos(Angle);
SourceArrayZr = SourceArrayZ + fCouchDistPerView*Angle + Zanchor;

DetectorArrayXr =  DetectorArrayX.*cos(Angle) + DetectorArrayY.*sin(Angle);
DetectorArrayYr = -DetectorArrayX.*sin(Angle) + DetectorArrayY.*cos(Angle);
DetectorArrayZr = DetectorArrayZ + fCouchDistPerView*Angle + Zanchor;
DetectorArrayXrArray = repmat(DetectorArrayXr,[SliceNum,1]);
DetectorArrayYrArray = repmat(DetectorArrayYr,[SliceNum,1]);
DetectorArrayZrArray = bsxfun(@plus, DetectorArrayZr, DetectorZArray);

h = plot3(parentH, [SourceArrayXr; DetectorArrayXrArray; SourceArrayXr],[SourceArrayYr; DetectorArrayYrArray; SourceArrayYr],[SourceArrayZr; DetectorArrayZrArray; SourceArrayZr], 'b.-');
% plot3([SourceArrayX; DetectorArrayX; DetectorArrayX; SourceArrayX],[SourceArrayY; DetectorArrayY; DetectorArrayY; SourceArrayY],[SourceArrayZ; DetectorArrayZTop; DetectorArrayZBot; SourceArrayZ]);

end
function DetecPlotRefresh()
global handles
delete(handles.Detect);
handles.Detect = DetecPlot(handles.ShowPanel, handles.GeoPara);
drawnow;
end
function AutoPlayTask(obj,event)
global handles
angle = handles.GeoPara.AngViewValue.Value + handles.GeoPara.AngViewStep.Value;
if (angle <= handles.GeoPara.AngViewMax.Value)
    handles.GeoPara.AngViewValue.Value = angle;
else
    handles.GeoPara.AngViewValue.Value = handles.GeoPara.AngViewMin.Value;
end
UIUpdateFromGeoPara('AngViewValue');
DetecPlotRefresh();
end
