function varargout = PreprocessImage(varargin)
% PREPROCESSIMAGE MATLAB code for PreprocessImage.fig
%      PREPROCESSIMAGE, by itself, creates a new PREPROCESSIMAGE or raises the existing
%      singleton*.
%
%      H = PREPROCESSIMAGE returns the handle to a new PREPROCESSIMAGE or the handle to
%      the existing singleton*.
%
%      PREPROCESSIMAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESSIMAGE.M with the given input arguments.
%
%      PREPROCESSIMAGE('Property','Value',...) creates a new PREPROCESSIMAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PreprocessImage_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PreprocessImage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PreprocessImage

% Last Modified by GUIDE v2.5 20-Jul-2016 16:29:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PreprocessImage_OpeningFcn, ...
                   'gui_OutputFcn',  @PreprocessImage_OutputFcn, ...
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


% --- Executes just before PreprocessImage is made visible.
function PreprocessImage_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PreprocessImage (see VARARGIN)

set(handles.FishViewerSlider,'Enable','off');
set(handles.FishNameEditText,'Enable','off','String','');
set(handles.SliceNumberEditText,'Enable','off','String','');
cla(handles.FishViewerAxes)
set(handles.FishViewerAxes,'Visible','off');
set(handles.UpButton,'Enable','off');
set(handles.DownButton,'Enable','off');
set(handles.LeftButton,'Enable','off');
set(handles.RightButton,'Enable','off');
set(handles.SaveButton,'Enable','off');
set(handles.RotatorLinkButton,'Enable','off');
set(handles.checkbox1,'Enable','off');
set(handles.AngleRetrieveButton,'Enable','off');
set(handles.RotationProcessButton,'Enable','off');
set(handles.AngleEditField,'Enable','off');

setappdata(0,'rotationAngle',[]);
handles.ViewOption = 0;
handles.allowFunction = 0;
handles.firstTime = 1;
% Choose default command line output for PreprocessImage
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PreprocessImage wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PreprocessImage_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ImageTypes = {...
    '*.dcm',...
    'DICOM Files (*.dcm)';...
    '*.jpg;*.png;*.tif',...
    'Image Files (*.jpg,*.png,*.tif)';...
    '*.img;*.hdr;*.analyze75',...
    'NIFTI Files (*.img,*.hdr,*.analyze75)'};
[FileName, PathName, FilterIndex] = uigetfile(ImageTypes,'Select Image File');
switch FilterIndex
    case 1
        I = squeeze(dicomread([PathName filesep FileName]));
    case 2
        I = load_nii([PathName filesep FileName]);
        I = I.img;
    case 3
        I = imread([PathName filesep FileName]);
end

if FilterIndex ~= 0
    
    set(handles.FishNameEditText,'String',FileName)
    im = mat2gray(I,[-2^15 2^15]);
    im = flip(permute(im, [1 3 2]),3);
    
    sip = squeeze(max(im,[],3));
    imshow(sip,'Parent',handles.FishViewerAxes,'Border','tight');
    colormap(bone)
    
    if handles.firstTime
        set(handles.UpButton,'Enable','on');
        set(handles.DownButton,'Enable','on');
        set(handles.LeftButton,'Enable','on');
        set(handles.RightButton,'Enable','on');
        set(handles.SaveButton,'Enable','on');
        set(handles.RotatorLinkButton,'Enable','on');
        set(handles.checkbox1,'Enable','on');
        set(handles.AngleRetrieveButton,'Enable','on');
        set(handles.RotationProcessButton,'Enable','on');
        set(handles.AngleEditField,'Enable','on');
        handles.firstTime = 0;
    end
    
    handles.FishImage = im;
    guidata(hObject, handles);
end

% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

save_folder = uigetdir;
if ischar(save_folder)
    allfiles = dir(save_folder);
    FishName = get(handles.FishNameEditText,'String');
    if sum(ismember(FishName,{allfiles.name})) > 0
        
        n=1;
        tempFishName = [FishName ' (' int2str(n) ')'];
        
        while sum(ismember(tempFishName,{allfiles.name})) > 0
            n = n+1;
            tempFishName(end-1) = int2str(n);
        end
        
        FishName = tempFishName;
    end
    
    [x,y,z] = size(handles.FishImage);
    dicomwrite(reshape(int16(handles.FishImage*2^16-2^15),[x y 1 z]),[save_folder filesep FishName]);
end

% --- Executes on button press in UpButton.
function UpButton_Callback(hObject, eventdata, handles)
% hObject    handle to UpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FishImage = flip(permute(handles.FishImage, [2 1 3]),2);

if handles.ViewOption == 0 
    
    sipDummy = max(handles.FishImage,[],3);
    sip = squeeze(sipDummy);
    imshow(sip,'Border','tight','Parent',handles.FishViewerAxes);
    colormap(bone);
    
else
    
    [~,~,thirdDim] = size(handles.FishImage);
    imshow(squeeze(handles.FishImage(:,:,round(thirdDim/2))));
    set(handles.SliceNumberEditText,'String',int2str(round(thirdDim/2)));

    set(handles.FishViewerSlider, 'SliderStep', [0.5/(round(thirdDim/2)-1) 1/(round(thirdDim/2)-1)*10], 'max', thirdDim, 'min', 1, 'Value',...
        round(thirdDim/2));
    
end

guidata(hObject, handles);
% --- Executes on button press in DownButton.
function DownButton_Callback(hObject, eventdata, handles)
% hObject    handle to DownButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FishImage = flip(permute(handles.FishImage, [2 1 3]),1);

if handles.ViewOption == 0 
    
    sipDummy = max(handles.FishImage,[],3);
    sip = squeeze(sipDummy);
    imshow(sip,'Border','tight','Parent',handles.FishViewerAxes);
    colormap(bone);
    
else
    
    [~,~,thirdDim] = size(handles.FishImage);
    imshow(squeeze(handles.FishImage(:,:,round(thirdDim/2))));
    set(handles.SliceNumberEditText,'String',int2str(round(thirdDim/2)));

    set(handles.FishViewerSlider, 'SliderStep', [0.5/(round(thirdDim/2)-1) 1/(round(thirdDim/2)-1)*10], 'max', thirdDim, 'min', 1, 'Value',...
        round(thirdDim/2));
    
end

guidata(hObject, handles);

% --- Executes on button press in LeftButton.
function LeftButton_Callback(hObject, eventdata, handles)
% hObject    handle to LeftButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FishImage = flip(permute(handles.FishImage, [1 3 2]),2);

if handles.ViewOption == 0 
    
    sipDummy = max(handles.FishImage,[],3);
    sip = squeeze(sipDummy);
    imshow(sip,'Border','tight','Parent',handles.FishViewerAxes);
    colormap(bone);
    
else
    
    [~,~,thirdDim] = size(handles.FishImage);
    imshow(squeeze(handles.FishImage(:,:,round(thirdDim/2))));
    set(handles.SliceNumberEditText,'String',int2str(round(thirdDim/2)));

    set(handles.FishViewerSlider, 'SliderStep', [0.5/(round(thirdDim/2)-1) 1/(round(thirdDim/2)-1)*10], 'max', thirdDim, 'min', 1, 'Value',...
        round(thirdDim/2));
    
end

guidata(hObject, handles);

% --- Executes on button press in RightButton.
function RightButton_Callback(hObject, eventdata, handles)
% hObject    handle to RightButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FishImage = flip(permute(handles.FishImage, [1 3 2]),3);

if handles.ViewOption == 0 
    
    sipDummy = max(handles.FishImage,[],3);
    sip = squeeze(sipDummy);
    imshow(sip,'Border','tight','Parent',handles.FishViewerAxes);
    colormap(bone);
    
else
    
    [~,~,thirdDim] = size(handles.FishImage);
    imshow(squeeze(handles.FishImage(:,:,round(thirdDim/2))));
    set(handles.SliceNumberEditText,'String',int2str(round(thirdDim/2)));

    set(handles.FishViewerSlider, 'SliderStep', [0.5/(round(thirdDim/2)-1) 1/(round(thirdDim/2)-1)*10], 'max', thirdDim, 'min', 1, 'Value',...
        round(thirdDim/2));
    
end

guidata(hObject, handles);

% --- Executes on slider movement.
function FishViewerSlider_Callback(hObject, eventdata, handles)
% hObject    handle to FishViewerSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currslice = round(get(hObject,'Value'));
imshow(squeeze(handles.FishImage(:,:,currslice)),'Parent',handles.FishViewerAxes);
set(handles.SliceNumberEditText,'String',int2str(currslice));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function FishViewerSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FishViewerSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ViewOption = get(hObject,'Value');
if handles.ViewOption
    
    [~,~,thirdDim] = size(handles.FishImage);
    
    set(handles.FishViewerSlider,'Enable','on');
    set(handles.FishViewerSlider, 'SliderStep', [0.5/(round(thirdDim/2)-1)...
        1/(round(thirdDim/2)-1)*10], 'max', thirdDim, 'min', 1, 'Value',...
        round(thirdDim/2));

    imshow(squeeze(handles.FishImage(:,:,round(thirdDim/2))),'Parent',handles.FishViewerAxes);
    set(handles.SliceNumberEditText,'String',int2str(round(thirdDim/2)));
    
else

    set(handles.FishViewerSlider,'Enable','off');
    set(handles.SliceNumberEditText,'String','');
    sip = squeeze(max(handles.FishImage,[],3));
    imshow(sip,'Parent',handles.FishViewerAxes,'Border','tight');
    colormap(bone);
    
end

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox1

function SliceNumberEditText_Callback(hObject, eventdata, handles)
% hObject    handle to SliceNumberEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SliceNumberEditText as text
%        str2double(get(hObject,'String')) returns contents of SliceNumberEditText as a double


% --- Executes during object creation, after setting all properties.
function SliceNumberEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliceNumberEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in RotatorLinkButton.
function RotatorLinkButton_Callback(hObject, eventdata, handles)
% hObject    handle to RotatorLinkButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageHandle = get(handles.FishViewerAxes,'Children');
setappdata(0,'fishSlice',imageHandle.CData);
Rotator

% --- Executes on button press in BackButton.
function BackButton_Callback(hObject, eventdata, handles)
% hObject    handle to BackButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(PreprocessImage)
close(Rotator)
FishCuT

function FishNameEditText_Callback(hObject, eventdata, handles)
% hObject    handle to FishNameEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FishNameEditText as text
%        str2double(get(hObject,'String')) returns contents of FishNameEditText as a double


% --- Executes during object creation, after setting all properties.
function FishNameEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FishNameEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AngleEditField_Callback(hObject, eventdata, handles)
% hObject    handle to AngleEditField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AngleEditField as text
%        str2double(get(hObject,'String')) returns contents of AngleEditField as a double


% --- Executes during object creation, after setting all properties.
function AngleEditField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AngleEditField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AngleRetrieveButton.
function AngleRetrieveButton_Callback(hObject, eventdata, handles)
% hObject    handle to AngleRetrieveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
angle = getappdata(0,'rotationAngle');
if ~isempty(angle)
    set(handles.AngleEditField,'String',num2str(angle));
end

% --- Executes on button press in RotationProcessButton.
function RotationProcessButton_Callback(hObject, eventdata, handles)
% hObject    handle to RotationProcessButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
angle = str2double(get(handles.AngleEditField,'String'));
fishim = handles.FishImage;
newfishim = imrotate(fishim,angle);
handles.FishImage = newfishim;

if handles.ViewOption == 0 
    
    sipDummy = max(handles.FishImage,[],3);
    sip = squeeze(sipDummy);
    imshow(sip,'Border','tight','Parent',handles.FishViewerAxes);
    colormap(bone);
    
else
    
    [~,~,thirdDim] = size(handles.FishImage);
    imshow(squeeze(handles.FishImage(:,:,round(thirdDim/2))));
    set(handles.SliceNumberEditText,'String',int2str(round(thirdDim/2)));

    set(handles.FishViewerSlider, 'SliderStep', [0.5/(round(thirdDim/2)-1) 1/(round(thirdDim/2)-1)*10], 'max', thirdDim, 'min', 1, 'Value',...
        round(thirdDim/2));
    
end
guidata(hObject, handles);
