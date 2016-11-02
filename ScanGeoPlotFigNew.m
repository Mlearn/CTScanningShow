function ScanGeoPlotFigNew()
% SCANGEOPLOTFIGNEW CT扫描结构示意图
% Created by An Mou
% Oct. 13, 2016

global handles

GeoParaCreate(); % 创建参数
BodyParaCreate();
UI_Build();      % 创建UI组件
GeoParaConectUI(); %连接参数和组件
UIConectGeoPara(); %连接组件和参数
UIValueInit();     %初始化UI组件中的值
TotalPlot();       %初步画图

handles.AutoPlayhandle = timer('Name','CircleTimer', 'TimerFcn', @AutoPlayTask, 'Period',0.2, 'ExecutionMode','fixedspacing');
handles.AutoPlayStatus = 'Original';
end
function UI_Build()
global handles
scrsz = get(0,'ScreenSize');

main_fig_size_position = [scrsz(1)+floor(0.3*scrsz(3)), scrsz(2)+floor(0.3*scrsz(4)), floor(0.5*scrsz(3)), floor(0.55*scrsz(4))];
handles.main_fig = figure('Tag','main_fig','Position',main_fig_size_position, 'Resize','on');
handles.ShowPanel = axes('visible','on','Position',[0 0 0.7 1],'XTick',[],'YTick',[],'ZTick',[],'Children',[],'Box','on','PlotBoxAspectRatio',[1,1,0.5], 'NextPlot','add');%,'CameraPosition',[1,0.1,0.05]);
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
% ==探测器层数
handles.UI.SliceNumName.handle=uicontrol(handles.main_fig,'Style','text','String','SliceNum:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(2,:));
handles.UI.SliceNumEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','32',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(2,:),'callback',{@Edit_GeometryReset_Callback, 'SliceNumEdit'});
% ==探测器一层厚度
handles.UI.SliceWidthName.handle=uicontrol(handles.main_fig,'Style','text','String','SliceWidth:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(3,:));
handles.UI.SliceWidthEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','0.625',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(3,:),'callback',{@Edit_GeometryReset_Callback, 'SliceWidthEdit'});
handles.UI.SliceWidthUnit.handle=uicontrol(handles.main_fig,'Style','text','String','mm',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(3,:)+[xW2 0 -xW2+0.06 0]);
% ==螺距
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
% ==投影角度范围下限
handles.UI.AngViewMinName.handle=uicontrol(handles.main_fig,'Style','text','String','ViewMin:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(6,:));
handles.UI.AngViewMinEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','0',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(6,:),'callback',{@Edit_GeometryReset_Callback, 'AngViewMinEdit'});
handles.UI.AngViewMinUnit.handle=uicontrol(handles.main_fig,'Style','text','String','π',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(6,:)+[xW2 0 -xW2+0.04 0]);
% ==投影角度范围上限
handles.UI.AngViewMaxName.handle=uicontrol(handles.main_fig,'Style','text','String','ViewMax:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(7,:));
handles.UI.AngViewMaxEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','5',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(7,:),'callback',{@Edit_GeometryReset_Callback, 'AngViewMaxEdit'});
handles.UI.AngViewMaxUnit.handle=uicontrol(handles.main_fig,'Style','text','String','π',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(7,:)+[xW2 0 -xW2+0.04 0]);
% ==投影角度步长
handles.UI.AngViewStepName.handle=uicontrol(handles.main_fig,'Style','text','String','ViewStep:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(8,:));
handles.UI.AngViewStepEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','0.01',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(8,:),'callback',{@Edit_GeometryReset_Callback, 'AngViewStepEdit'});
handles.UI.AngViewStepUnit.handle=uicontrol(handles.main_fig,'Style','text','String','π',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(8,:)+[xW2 0 -xW2+0.04 0]);
% ==投影角度
handles.UI.AngViewName.handle=uicontrol(handles.main_fig,'Style','text','String','ProjView:',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoName(9,:));
handles.UI.AngViewEdit.handle=uicontrol(handles.main_fig,'Style','edit','String','0',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(9,:),'callback',{@Scan_Callback, 'AngViewEdit'});
handles.UI.AngViewUnit.handle=uicontrol(handles.main_fig,'Style','text','String','π',...
    'HorizontalAlignment','left','FontWeight','normal','FontUnits','normalized','FontSize',0.7,'Units','normalized',...
    'Position',PoValu(9,:)+[xW2 0 -xW2+0.04 0]);

% ==
handles.UI.AutoPlayButton.Handle = uicontrol('Style', 'pushbutton', 'String', 'Play', 'Units','normalized',...
        'Position', [0.72 0.04 0.08 0.06], 'Callback', {@Button_AutoPlay_Callback, 'Play'});
handles.UI.AutoPlayButton.Handle = uicontrol('Style', 'pushbutton', 'String', 'Pause', 'Units','normalized',...
        'Position', [0.8 0.04 0.08 0.06], 'Callback', {@Button_AutoPlay_Callback, 'Pause'});
handles.UI.AutoPlayButton.Handle = uicontrol('Style', 'pushbutton', 'String', 'Reset', 'Units','normalized',...
        'Position', [0.88 0.04 0.08 0.06], 'Callback', {@Button_AutoPlay_Callback, 'Reset'});
% ==加载图像
handles.UI.AutoPlayButton.Handle = uicontrol('Style', 'pushbutton', 'String', 'Load', 'Units','normalized',...
        'Position', [0.88 0.24 0.08 0.06], 'Callback', @Loadfiles_Callback);
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
handles.GeoPara.ChannelNumPar.Value = floor(1665/Samp);
handles.GeoPara.MidChannelPar.Value = 8.265021361983943e+02/Samp;
handles.GeoPara.ChannelParSpace.Value = 0.305497508156251*Samp;
end
function BodyParaCreate()
global handles
handles.BodyPara.FOV.Value = 300;
handles.BodyPara.XCenter.Value = 0;
handles.BodyPara.YCenter.Value = 0;
handles.BodyPara.XPixelNum = 512;
handles.BodyPara.YPixelNum = 512;
handles.BodyPara.SliceThick.Value = 2;
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
% 每一次有一个UI变动，就需要调用一次该函数
global handles
ParaName = [];
if isfield(handles.UI.(UIName),'relatePara')
    ParaName = handles.UI.(UIName).relatePara{1};
    if strcmp('edit', handles.UI.(UIName).handle.Style) % 判断UI为可编辑输入框
        GeoParaValueSet(ParaName, str2double(get(handles.UI.(UIName).handle,'string')));
    end
    if strcmp('slider', handles.UI.(UIName).handle.Style) % 判断UI为滑块
        GeoParaValueSet(ParaName, get(handles.UI.(UIName).handle,'Value'));
    end
end
end
function UIUpdateFromGeoPara(ParaName)
% 改动了GeoPara之后需要调用该函数一次
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
    GeoParaCreate(); % 创建参数
    UIValueInit();     %初始化UI组件中的值
    TotalPlot();       %初步画图
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
CylinderPlot(handles.ShowPanel, GeoParameter);
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
axes(ParentH);%若想累加画图效果有效，必须先指定axes()，而不能直接在plot3函数中设定父axes
h = plot3(Xarray,Yarray,Zarray,'r-');

end
function h = PlanePlot(ParentH, Zposi, GeoParameter)

axes(ParentH);% hold on
R = GeoParameter.SFOV.Value/2;
ImgZ = Zposi;
%一种画透明平面方法（只能是四方的吗？）
% xy = [-hb, hb];
% surf(xy,xy, ones(length(xy),length(xy))*ImgZ); % 画出重建平面,不加Parent时候才能全部都画出来
% alpha(0.5);shading interp;
% 另一种画圆，只能在坐标z=0处
% h = rectangle('Position',[-hb,-hb,2*hb,2*hb],'Curvature',[1,1],  'FaceColor','r'); %只能做在z=0的平面
%另一种画平面的方法，可以是任意多边形
% vv = [-R,-R,ImgZ; R,-R,ImgZ;R,R,ImgZ;-R,R,ImgZ];
% patch('vertices', vv, 'faces', [1,2,3,4], 'FaceColor', 'y','FaceAlpha',0.5);
samp = 100;
x = R*sin(linspace(0,2*pi,samp));
y = R*cos(linspace(0,2*pi,samp));
zc = ImgZ*ones(1, samp);
vv = [x',y',zc'];
patch('vertices', vv, 'faces', 1:samp, 'FaceColor', 'r','FaceAlpha',0.9, 'Linestyle','-');%中间FOV
end
function h = ImagesPlot(ParentH, GeoParameter)
global handles
axes(ParentH);% hold on
[xNum,yNum,zNum] = size(handles.Data.Images);
R = handles.BodyPara.FOV.Value;
XCenter = handles.BodyPara.XCenter.Value;
YCenter = handles.BodyPara.YCenter.Value;
% XPixelNum = handles.BodyPara.XPixelNum.Value;
% YPixelNum = handles.BodyPara.YPixelNum.Value;
SliceThick = handles.BodyPara.SliceThick.Value;
x = XCenter+linspace(-R,R,yNum);
y = YCenter+linspace(-R,R,xNum);
z = SliceThick*(-(zNum-1)/2:(zNum-1)/2);
[x,y,z] = meshgrid(x,y,z);
zslice = 1:zNum;
Images = double(handles.Data.Images);
Images(Images>100) = 100;
Images(Images<-100) = -100;
Images = 100-Images;
Images(Images<-100) = NaN;
h = slice(x,y,z,Images,[],[],zslice);
shading interp;
colormap hot;
alpha(h, 0.3);
end
function h = CylinderPlot(ParentH, GeoParameter)
axes(ParentH);
R = GeoParameter.SFOV.Value/2;
SliceWidth = GeoParameter.SliceWidth.Value;
SliceNum = GeoParameter.SliceNum.Value;
DetectorZ = SliceWidth*SliceNum;
Pitch = GeoParameter.Pitch.Value;
fCouchDistPerView = DetectorZ*Pitch/2/pi;
AngViewMin = GeoParameter.AngViewMin.Value;
AngViewMax = GeoParameter.AngViewMax.Value;
samp = 100;
x = R*sin(linspace(0,2*pi,samp));
y = R*cos(linspace(0,2*pi,samp));
zb = pi*AngViewMin*fCouchDistPerView*ones(1, samp);
zt = pi*AngViewMax*fCouchDistPerView*ones(1, samp);
vv = [0,0,zb(1); x',y',zb'; 0,0,zt(1); x',y',zt'];
pOder = [ones(1, samp-1)',1+(1:samp-1)',2+(1:samp-1)',ones(1, samp-1)';
         bsxfun(@plus, (0:samp-2)', [2,3,samp+4,samp+3]);
         (samp+2)*ones(1, samp-1)',samp+2+(1:samp-1)',samp+3+(1:samp-1)',(samp+2)*ones(1, samp-1)'];
h = patch('vertices', vv, 'faces', pOder, 'FaceColor', 'y','FaceAlpha',0.5, 'Linestyle','none');%圆柱侧边
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

% axes(parentH);%添加该语句之后，图像总是刷新，在windows预览窗口下会频繁闪烁
vv = [SourceArrayXr',SourceArrayYr',SourceArrayZr';DetectorArrayXr',DetectorArrayYr',DetectorArrayZr'+DetectorZArray(1);DetectorArrayXr',DetectorArrayYr',DetectorArrayZr'+DetectorZArray(end)];
pOder = [1:ChannelNumPar,ChannelNumPar+(ChannelNumPar:-1:1);
         1:ChannelNumPar,2*ChannelNumPar+(ChannelNumPar:-1:1);
         ChannelNumPar+(1:ChannelNumPar),2*ChannelNumPar+(ChannelNumPar:-1:1)];
h = patch('Parent',parentH,'vertices', vv, 'faces', pOder, 'FaceColor', 'g','FaceAlpha',0.5, 'Linestyle','--');
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

function Loadfiles_Callback(hObject, callbackdata)
global handles
[FileNames,PathName] = getFileNames([]);
if PathName(1)==0
    return;
end
for i=2:length(FileNames)
    if (~strcmp(FileNames{i}(end),FileNames{1}(end)))&&~isempty(strfind(FileNames{i},'.'))
        errordlg('Different Image Format!');
        return;
    end
end

    FolderName = PathName;

FigureFullNames=fullfile(FolderName,FileNames);
if isempty(find(FileNames{1}=='.', 1))
    switch questdlg('No extern name!Open as dicom files?','warn','yes','no','yes')
        case 'yes'
            [CurrentValue.Imgdata, CurrentValue.ImageTotalNumber] = loadimages(FigureFullNames,'dcm');
        case 'no'
            return;
    end
else
    switch FileNames{1}(end-2:end)
        case 'svn'
            for im=1:length(FileNames)
                filename=fullfile(PathName,FileNames(im));
                filename=char(filename);
                fp=fopen(filename,'r');
                m=fread(fp,1,'int');
                n=fread(fp,1,'int');
                minclim=fread(fp,1,'int');
                maxclim=fread(fp,1,'int');
                CurrentValue.Imgdata(:,:,im)=fread(fp,[m n],'double');
                %#function imtool
                imtool(CurrentValue.Imgdata(:,:,im),[minclim maxclim]);
            end
        case {'mat','MAT'}
            [CurrentValue.Imgdata, CurrentValue.ImageTotalNumber] = loadimages(FigureFullNames,'mat');
        otherwise %'dcm'
            [CurrentValue.Imgdata, CurrentValue.ImageTotalNumber] = loadimages(FigureFullNames,'dcm');
    end
end
GeoParameter = handles.GeoPara;
ImagesPlot(handles.ShowPanel, GeoParameter);
end
function [FileNames,PathName] = getFileNames(foldername)

if isempty(foldername) && ~exist('newpath','var')
    foldername=cd;
    foldername = strcat(foldername,'\');
end
[FileNames,PathName]=uigetfile({'*.*','All Files';'*.mat','MAT-files(*.mat)';'*.dcm','DICOM Files(*.dcm)'},'Load Mat/DCM-Type images', foldername,'MultiSelect','on');
tf = iscell(FileNames);
if tf ~=1
%     if FileNames(1) == 0 || FileNames(2) == 0 
    if isequal(FileNames,0)
         return;
    end
    FileNames = {FileNames}; 
end
end

function [Imgdata, ImageNum] = loadimages(FileNames,Format)

global handles
ImageNum = length(FileNames);
if strcmp(Format,'mat')
    matstruct=load(FileNames{1});
    temp=double(matstruct.img);
    [iHei, iWid] = size(temp);
    Imgdata = zeros(iHei,iWid,ImageNum);
    Imgdata(:,:,1) = temp;
    h = waitbar(1/ImageNum,'Please wait...');
    for i = 2:ImageNum
        matstruct=load(FileNames{i});
        Imgdata(:,:,i)=double(matstruct.img);
        waitbar(i / ImageNum);
    end
    close(h);
elseif strcmp(Format,'dcm')
    info=dicominfo(FileNames{1});
    temp=dicomread(FileNames{1});
    [iHei, iWid] = size(temp);
    Imgdata = zeros(iHei,iWid, ImageNum);
    temp = single(temp)*info.RescaleSlope+info.RescaleIntercept;
    Imgdata(:,:,1) = temp;
    h = waitbar(1/ImageNum,'Please wait...');
    for i = 2:ImageNum
        infotemp=dicominfo(FileNames{i});
        Imgdata(:,:,i) = dicomread(FileNames{i});
        Imgdata(:,:,i)= single(Imgdata(:,:,i))*infotemp.RescaleSlope+infotemp.RescaleIntercept;
        waitbar(i / ImageNum);
    end
    close(h);
end
handles.Data.Images = Imgdata;

end