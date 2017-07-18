function varargout = ChangeParameters(varargin)
% CHANGEPARAMETERS MATLAB code for ChangeParameters.fig
%      CHANGEPARAMETERS, by itself, creates a new CHANGEPARAMETERS or raises the existing
%      singleton*.
%
%      H = CHANGEPARAMETERS returns the handle to a new CHANGEPARAMETERS or the handle to
%      the existing singleton*.
%
%      CHANGEPARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHANGEPARAMETERS.M with the given input arguments.
%
%      CHANGEPARAMETERS('Property','Value',...) creates a new CHANGEPARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChangeParameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChangeParameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ChangeParameters

% Last Modified by GUIDE v2.5 28-Jun-2017 15:32:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChangeParameters_OpeningFcn, ...
                   'gui_OutputFcn',  @ChangeParameters_OutputFcn, ...
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


% --- Executes just before ChangeParameters is made visible.
function ChangeParameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChangeParameters (see VARARGIN)

[PATHSTR,~,~] = fileparts(mfilename('fullpath'));
fileID = fopen([PATHSTR filesep 'ImagingParameters.txt']);
C = textscan(fileID,'%s %f');
fclose(fileID);

Params = strtrim(cellstr(num2str(C{2})));

set(handles.ResolutionEdit,'String',Params{1});
set(handles.SlopeEdit,'String',Params{2});
set(handles.InterceptEdit,'String',Params{3});
set(handles.ThresholdMultiplierEdit,'String',Params{4});
handles.Params = Params;

% Choose default command line output for ChangeParameters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ChangeParameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ChangeParameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ResolutionEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ResolutionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ResolutionEdit as text
%        str2double(get(hObject,'String')) returns contents of ResolutionEdit as a double


% --- Executes during object creation, after setting all properties.
function ResolutionEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ResolutionEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SlopeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SlopeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SlopeEdit as text
%        str2double(get(hObject,'String')) returns contents of SlopeEdit as a double


% --- Executes during object creation, after setting all properties.
function SlopeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlopeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InterceptEdit_Callback(hObject, eventdata, handles)
% hObject    handle to InterceptEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InterceptEdit as text
%        str2double(get(hObject,'String')) returns contents of InterceptEdit as a double


% --- Executes during object creation, after setting all properties.
function InterceptEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InterceptEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EndButton.
function EndButton_Callback(hObject, eventdata, handles)
% hObject    handle to EndButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Param1 = get(handles.ResolutionEdit, 'String');
Param2 = get(handles.SlopeEdit, 'String');
Param3 = get(handles.InterceptEdit, 'String');
Param4 = get(handles.ThresholdMultiplierEdit, 'String');

[PATHSTR,~,~] = fileparts(mfilename('fullpath'));
fileID = fopen([PATHSTR filesep 'ImagingParameters.txt'],'w');
C = {['Resolution: ' Param1];['Slope: ' Param2];['Intercept: ' Param3]; ['ThresholdMultiplier: ' Param4]};
[nrows, ~] = size(C);
for row = 1:nrows
    temp_str = C{row,:};
    fprintf(fileID, '%s\n', temp_str);
end

fclose(fileID);
close(ChangeParameters);


% --- Executes on button press in UndoButton.
function UndoButton_Callback(hObject, eventdata, handles)
% hObject    handle to UndoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ResolutionEdit,'String',handles.Params{1});
set(handles.SlopeEdit,'String',handles.Params{2});
set(handles.InterceptEdit,'String',handles.Params{3});
set(handles.ThresholdMultiplierEdit,'String',handles.Params{4});

guidata(hObject, handles);



function ThresholdMultiplierEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ThresholdMultiplierEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThresholdMultiplierEdit as text
%        str2double(get(hObject,'String')) returns contents of ThresholdMultiplierEdit as a double


% --- Executes during object creation, after setting all properties.
function ThresholdMultiplierEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThresholdMultiplierEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
