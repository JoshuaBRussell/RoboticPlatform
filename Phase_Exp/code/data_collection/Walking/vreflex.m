function varargout = vreflex(varargin)
% VREFLEX MATLAB code for vreflex.fig
%      VREFLEX, by itself, creates a new VREFLEX or raises the existing
%      singleton*.
%
%      H = VREFLEX returns the handle to a new VREFLEX or the handle to
%      the existing singleton*.
%
%      VREFLEX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VREFLEX.M with the given input arguments.
%
%      VREFLEX('Property','Value',...) creates a new VREFLEX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vreflex_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vreflex_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vreflex

% Last Modified by GUIDE v2.5 21-Aug-2019 19:57:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vreflex_OpeningFcn, ...
                   'gui_OutputFcn',  @vreflex_OutputFcn, ...
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


% --- Executes just before vreflex is made visible.
function vreflex_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vreflex (see VARARGIN)

% Choose default command line output for vreflex
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vreflex wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vreflex_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setup.
function setup_Callback(hObject, eventdata, handles)
% hObject    handle to setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amplitude=0;
assignin('base','amplitude',amplitude);
phi=3.14;
assignin('base','phi',phi);
off_TA = 0;
assignin('base','off_TA',off_TA);
off_SOL = 0;
assignin('base','off_SOL',off_SOL);
TA_MVC = 1;
assignin('base','TA_MVC',TA_MVC);
SOL_MVC = 1;
assignin('base','SOL_MVC',SOL_MVC);
TMVC = str2double(get(handles.act,'String'))/100;
assignin('base','TMVC',TMVC);
Ts=0.0005;
assignin('base','Ts',Ts);
trial=0;
assignin('base','trial',trial);
range=0;
assignin('base','range',range);
DPbase=0;
assignin('base','DPbase',DPbase);
IEbase=0;
assignin('base','IEbase',IEbase);
forrange=5;
assignin('base','forrange',forrange);
maxtorquediff=50;
assignin('base','maxtorquediff',maxtorquediff);
f1_off=0;
f2_off=0;
f3_off=0;
f4_off=0;
f5_off=0;
f6_off=0;
assignin('base','f1_off',f1_off);
assignin('base','f2_off',f2_off);
assignin('base','f3_off',f3_off);
assignin('base','f4_off',f4_off);
assignin('base','f5_off',f5_off);
assignin('base','f6_off',f6_off);

% --- Executes on button press in update_mvc.
function update_mvc_Callback(hObject, eventdata, handles)
% hObject    handle to update_mvc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename,Pathname,Filtercopy]=uigetfile ;
% Now you can set the "string" property of the static textbox. 
if Filename~=0
    s=strcat(Pathname,Filename);
    name=get(handles.name,'String');
    r=strcat(name,'mvc_info.dat');
    copyfile(s,r)
    assignin('base','Filename',r);
    mvc_reflex
end

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
act = str2double(get(handles.act,'String'))/100;
assignin('base','TMVC',act);
trial = str2double(get(handles.tryno,'String'));
assignin('base','trial',trial);
if act ==0
    range=1;
else
    range=0.025
end 
assignin('base','range',range);
CoP=str2double(get(handles.copno,'String'));
weight=evalin('base','weight');
assignin('base','CoP',CoP);
 DP_cop = evalin('base','DP_cop');
 CoPbase=evalin('base','CoPbase');
 IE_cop = evalin('base','IE_cop');
 weight=evalin('base','weight');
 tweight=str2double(get(handles.weightpc,'String'))*0.01*weight;
 if trial==1
    
     IEF=IE_cop*tweight*0.01;
     DPF=(CoP+DP_cop)*(tweight)*0.01;
    
 else
   
     IEF=(CoP+IE_cop)*tweight*0.01;
     DPF=(DP_cop)*(tweight)*0.01;
     
 end
 
 assignin('base','DPF',DPF);
 assignin('base','IEF',IEF);
 assignin('base','tweight',tweight);
 
 

function act_Callback(hObject, eventdata, handles)
% hObject    handle to act (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of act as text
%        str2double(get(hObject,'String')) returns contents of act as a double


% --- Executes during object creation, after setting all properties.
function act_CreateFcn(hObject, eventdata, handles)
% hObject    handle to act (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tryno_Callback(hObject, eventdata, handles)
% hObject    handle to tryno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tryno as text
%        str2double(get(hObject,'String')) returns contents of tryno as a double


% --- Executes during object creation, after setting all properties.
function tryno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tryno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function name_Callback(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name as text
%        str2double(get(hObject,'String')) returns contents of name as a double


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename,Pathname,Filtercopy]=uigetfile ;
% Now you can set the "string" property of the static textbox. 
if Filename~=0
    s=strcat(Pathname,Filename);
    name=get(handles.name,'String');
    act=get(handles.act,'String');
    tryno=get(handles.tryno,'String');
    r=strcat(name,'_',act,'_',tryno,'.dat');
    copyfile(s,r)
end



function copno_Callback(hObject, eventdata, handles)
% hObject    handle to copno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of copno as text
%        str2double(get(hObject,'String')) returns contents of copno as a double


% --- Executes during object creation, after setting all properties.
function copno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to copno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function weightpc_Callback(hObject, eventdata, handles)
% hObject    handle to weightpc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of weightpc as text
%        str2double(get(hObject,'String')) returns contents of weightpc as a double


% --- Executes during object creation, after setting all properties.
function weightpc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to weightpc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
