function varargout = FishCuT(varargin)
% FISHCUT MATLAB code for FishCuT.fig
%      FISHCUT, by itself, creates a new FISHCUT or raises the existing
%      singleton*.
%
%      H = FISHCUT returns the handle to a new FISHCUT or the handle to
%      the existing singleton*.
%
%      FISHCUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FISHCUT.M with the given input arguments.
%
%      FISHCUT('Property','Value',...) creates a new FISHCUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FishCuT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FishCuT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FishCuT

% Last Modified by GUIDE v2.5 18-Jul-2017 11:52:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FishCuT_OpeningFcn, ...
                   'gui_OutputFcn',  @FishCuT_OutputFcn, ...
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


% --- Executes just before FishCuT is made visible.
function FishCuT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FishCuT (see VARARGIN)
thisPath = fileparts(mfilename('fullpath'));
backgroundImage = imread([thisPath filesep 'StartBackground.tif']);
ha = axes('units','normalized','position',[0 0 1 1]);
imshow(backgroundImage(:,:,1:3))
% Choose default command line output for FishCuT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FishCuT wait for user response (see UIRESUME)
% uiwait(handles.handle1);


% --- Outputs from this function are returned to the command line.
function varargout = FishCuT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in StartButton.
function StartButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles    structure with handles and user data (see GUIDATA)
ImageTypes = {...
    '*.dcm',...
    'DICOM Files (*.dcm)';...
    '*.jpg;*.png;*.tif',...
    'Image Files (*.jpg,*.png,*.tif)'};
 
[fileName, filePath, FilterIndex] = uigetfile(ImageTypes,'Select Image File');

parambool = ParameterCheck;

if ischar(fileName) && strcmp(parambool, 'OK')
    FishCuTScriptv1_1(filePath,fileName,FilterIndex);
elseif strcmp(parambool, 'Back')
    ChangeParameters
end

% --- Executes on button press in ExitButton.
function ExitButton_Callback(hObject, eventdata, handles)
% hObject    handle to ExitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(Rotator)
close(PreprocessImage)
close(FishCuT)

% --- Executes on button press in PreProcessButton.
function PreProcessButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreProcessButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PreprocessImage


% --- Executes on button press in handle1.
function ChangeImaging_Callback(hObject, eventdata, handles)
% hObject    handle to handle1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ChangeParameters
