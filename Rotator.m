function varargout = Rotator(varargin)
% ROTATOR MATLAB code for Rotator.fig
%      ROTATOR, by itself, creates a new ROTATOR or raises the existing
%      singleton*.
%
%      H = ROTATOR returns the handle to a new ROTATOR or the handle to
%      the existing singleton*.
%
%      ROTATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROTATOR.M with the given input arguments.
%
%      ROTATOR('Property','Value',...) creates a new ROTATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Rotator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Rotator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Rotator

% Last Modified by GUIDE v2.5 20-Jul-2016 15:52:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Rotator_OpeningFcn, ...
                   'gui_OutputFcn',  @Rotator_OutputFcn, ...
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


% --- Executes just before Rotator is made visible.
function Rotator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Rotator (see VARARGIN)
fishSlice = getappdata(0, 'fishSlice');
imshow(fishSlice,'Parent',handles.RotateSliceAxes,'Border','Tight')
handles.ofishSlice = fishSlice;
handles.fishSlice = fishSlice;

% Choose default command line output for Rotator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Rotator wait for user response (see UIRESUME)
% uiwait(handles.Rotator);


% --- Outputs from this function are returned to the command line.
function varargout = Rotator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ThresholdButton.
function ThresholdButton_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mijirun = 1;
classpath = javaclasspath('-all');
i=1;
while i <= size(classpath,1)
    [~,name,~] = fileparts(classpath{i});
    if strcmp(name,'z_spacing-0.0.1-SNAPSHOT')
        mijirun = 0;
        break
    end
    i = i+1;
end

if mijirun
    Miji(false);
end

Algorithm = 'Default';
fishSlice = handles.fishSlice;
fishSlice(fishSlice==0.5) = 0;
MIJ.createImage('SIP',uint16(fishSlice*2^16),true);
MIJ.run('Auto Threshold',['method=' Algorithm ' ignore_black ignore_white white']);
threshindxs = MIJ.getCurrentImage;
MIJ.run('Close All')
threshindxs(threshindxs~=0) = 1;
fishSlice(~logical(threshindxs)) = 0;
imshow(fishSlice,'Parent',handles.RotateSliceAxes,'Border','Tight')
handles.fishSlice = fishSlice;
guidata(hObject, handles);

% --- Executes on button press in BinarizeButton.
function BinarizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to BinarizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fishSlice = handles.fishSlice;
fishSlice(fishSlice==0.5) = 0;
fishSlice(fishSlice~=0) = 1;
fishSlice = medfilt2(fishSlice);
imshow(fishSlice,'Parent',handles.RotateSliceAxes,'Border','Tight')
handles.fishSlice = fishSlice;
guidata(hObject, handles);

% --- Executes on button press in FilterButton.
function FilterButton_Callback(hObject, eventdata, handles)
% hObject    handle to FilterButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fishSlice = handles.fishSlice;
connProp = bwconncomp(fishSlice,8);
pixelGroups = connProp.PixelIdxList;
pixelNums = cellfun(@numel,pixelGroups);
[~,maxGroup] = max(pixelNums);
maxIDs = pixelGroups{maxGroup};
imgIndxs = ones(size(fishSlice,1),size(fishSlice,2));
imgIndxs(maxIDs) = 0;
fishSlice(logical(imgIndxs)) = 0;
imshow(fishSlice,'Parent',handles.RotateSliceAxes,'Border','Tight')
handles.fishSlice = fishSlice;
guidata(hObject, handles);

% --- Executes on button press in BackButton.
function BackButton_Callback(hObject, eventdata, handles)
% hObject    handle to BackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(0,'rotationAngle',str2double(get(handles.RotateAngleTextField,'String')))
close(Rotator)

% --- Executes on button press in Compute_AngleButton.
function Compute_AngleButton_Callback(hObject, eventdata, handles)
% hObject    handle to Compute_AngleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fishSlice = handles.fishSlice;

if length(unique(fishSlice)) == 2
    fishSliceProps = regionprops(fishSlice, 'Orientation', 'MajorAxisLength', ...
    'MinorAxisLength', 'Eccentricity', 'Centroid');

    angle = fishSliceProps.Orientation;
    if angle > 45
        angle = 90-angle;
    elseif angle < -45
        angle = -90-angle;
    end
    
    setappdata(0,'rotationAngle',angle);
    set(handles.RotateAngleTextField,'String',num2str(angle));
    
    phi = linspace(0,2*pi,50);
    cosphi = cos(phi);
    sinphi = sin(phi);

    xbar = fishSliceProps.Centroid(1);
    ybar = fishSliceProps.Centroid(2);

    a = fishSliceProps.MajorAxisLength/2;
    b = fishSliceProps.MinorAxisLength/2;

    theta = pi*fishSliceProps.Orientation/180;
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];

    xy = [a*cosphi; b*sinphi];
    xy = R*xy;
    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;
    hold on
    axes(handles.RotateSliceAxes) 
    plot(x,y,'r','LineWidth',2);
end

% --- Executes on button press in RotateSliceButon.
function RotateSliceButon_Callback(hObject, eventdata, handles)
% hObject    handle to RotateSliceButon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ofishSlice = handles.ofishSlice;
fishSlice = handles.fishSlice;
angle = str2double(get(handles.RotateAngleTextField,'String'));

handles.ofishSlice = imrotate(ofishSlice,angle,'crop','bilinear');
fishSlice = imrotate(fishSlice,angle,'crop','bilinear');
handles.fishSlice = fishSlice;

imshow(fishSlice,'Parent',handles.RotateSliceAxes,'Border','Tight')
handles.fishSlice = fishSlice;
guidata(hObject, handles);

function RotateAngleTextField_Callback(hObject, eventdata, handles)
% hObject    handle to RotateAngleTextField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RotateAngleTextField as text
%        str2double(get(hObject,'String')) returns contents of RotateAngleTextField as a double


% --- Executes during object creation, after setting all properties.
function RotateAngleTextField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RotateAngleTextField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RevertButton.
function RevertButton_Callback(hObject, eventdata, handles)
% hObject    handle to RevertButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fishSlice = handles.ofishSlice;
handles.fishSlice = fishSlice;
imshow(fishSlice,'Parent',handles.RotateSliceAxes,'Border','Tight')
handles.fishSlice = fishSlice;
guidata(hObject, handles);
