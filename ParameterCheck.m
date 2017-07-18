function varargout = ParameterCheck(varargin)
% PARAMETERCHECK MATLAB code for ParameterCheck.fig
%      PARAMETERCHECK by itself, creates a new PARAMETERCHECK or raises the
%      existing singleton*.
%
%      H = PARAMETERCHECK returns the handle to a new PARAMETERCHECK or the handle to
%      the existing singleton*.
%
%      PARAMETERCHECK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETERCHECK.M with the given input arguments.
%
%      PARAMETERCHECK('Property','Value',...) creates a new PARAMETERCHECK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ParameterCheck_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ParameterCheck_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ParameterCheck

% Last Modified by GUIDE v2.5 28-Jun-2017 15:45:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ParameterCheck_OpeningFcn, ...
                   'gui_OutputFcn',  @ParameterCheck_OutputFcn, ...
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

% --- Executes just before ParameterCheck is made visible.
function ParameterCheck_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ParameterCheck (see VARARGIN)

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

% Choose default command line output for ParameterCheck
handles.output = 'Yes';

% Update handles structure
guidata(hObject, handles);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.
if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
         case 'title'
          set(hObject, 'Name', varargin{index+1});
         case 'string'
          set(handles.StaticTag, 'String', varargin{index+1});
        end
    end
end

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes ParameterCheck wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ParameterCheck_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in OKbutton.
function OKbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OKbutton (see GCBO)
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

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);

% --- Executes on button press in backButton.
function backButton_Callback(hObject, eventdata, handles)
% hObject    handle to backButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 'No';
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.figure1);
end    
    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.figure1);
end    



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
