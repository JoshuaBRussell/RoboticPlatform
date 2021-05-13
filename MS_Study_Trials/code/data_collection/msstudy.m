function varargout = msstudy(varargin)
% MSSTUDY MATLAB code for msstudy.fig
%      MSSTUDY, by itself, creates a new MSSTUDY or raises the existing
%      singleton*.
%
%      H = MSSTUDY returns the handle to a new MSSTUDY or the handle to
%      the existing singleton*.
%
%      MSSTUDY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSSTUDY.M with the given input arguments.
%
%      MSSTUDY('Property','Value',...) creates a new MSSTUDY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before msstudy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to msstudy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help msstudy

% Last Modified by GUIDE v2.5 31-Jan-2020 19:38:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @msstudy_OpeningFcn, ...
    'gui_OutputFcn',  @msstudy_OutputFcn, ...
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


% --- Executes just before msstudy is made visible.
function msstudy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to msstudy (see VARARGIN)

% Choose default command line output for msstudy
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes msstudy wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = msstudy_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function kdp_Callback(hObject, eventdata, handles)
% hObject    handle to kdp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kdp as text
%        str2double(get(hObject,'String')) returns contents of kdp as a double


% --- Executes during object creation, after setting all properties.
function kdp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kdp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kie_Callback(hObject, eventdata, handles)
% hObject    handle to kie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kie as text
%        str2double(get(hObject,'String')) returns contents of kie as a double


% --- Executes during object creation, after setting all properties.
function kie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DPStanding.
function DPStanding_Callback(hObject, eventdata, handles)
% hObject    handle to DPStanding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DPStanding


% --- Executes on button press in IEstanding.
function IEstanding_Callback(hObject, eventdata, handles)
% hObject    handle to IEstanding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IEstanding


% --- Executes on button press in DPwalking.
function DPwalking_Callback(hObject, eventdata, handles)
% hObject    handle to DPwalking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DPwalking


% --- Executes on button press in IEwalking.
function IEwalking_Callback(hObject, eventdata, handles)
% hObject    handle to IEwalking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of IEwalking


% --- Executes on button press in Set.
function Set_Callback(hObject, eventdata, handles)
% hObject    handle to Set (see GCBO)
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
TMVC = 0;
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
f7_off=0;
f8_off=0;
kdp=1000;
kie=1000;
assignin('base','kdp',kdp);
assignin('base','kie',kie);
assignin('base','f1_off',f1_off);
assignin('base','f2_off',f2_off);
assignin('base','f3_off',f3_off);
assignin('base','f4_off',f4_off);
assignin('base','f5_off',f5_off);
assignin('base','f6_off',f6_off);
assignin('base','f7_off',f7_off);
assignin('base','f8_off',f8_off);

f1_off_L=0;
f2_off_L=0;
f3_off_L=0;
f4_off_L=0;
f5_off_L=0;
f6_off_L=0;
f7_off_L=0;
f8_off_L=0;
assignin('base','f1_off_L',f1_off_L);
assignin('base','f2_off_L',f2_off_L);
assignin('base','f3_off_L',f3_off_L);
assignin('base','f4_off_L',f4_off_L);
assignin('base','f5_off_L',f5_off_L);
assignin('base','f6_off_L',f6_off_L);
assignin('base','f7_off_L',f7_off_L);
assignin('base','f8_off_L',f8_off_L);

% --- Executes on button press in updatemvc.
function updatemvc_Callback(hObject, eventdata, handles)
% hObject    handle to updatemvc (see GCBO)
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
if(get(handles.cop,'Value')==1)
    kdp = str2double(get(handles.kdp,'String'));
    kie=str2double(get(handles.kie,'String'));
    
    assignin('base','kdp',kdp);
    assignin('base','kie',kie);
    if(get(handles.visual,'Value')==1)
        visual=0;
        assignin('base','visual',visual);
    elseif(get(handles.visual,'Value')==0)
                visual=1;
        assignin('base','visual',visual);
    end
    ttime=2000*str2double(get(handles.ttime,'String'));
    assignin('base','ttime',ttime);
end
if(get(handles.standing,'Value')==1)
    if(get(handles.dpstanding,'Value')==1)
        act = 0;
        assignin('base','TMVC',act);
        trial = 1;
        assignin('base','trial',trial);
        if act ==0
            range=1;
        else
            range=0.025
        end
        assignin('base','range',range);
        CoP=0;
        weight=evalin('base','weight');
        assignin('base','CoP',CoP);
        DP_cop = evalin('base','DP_cop');
        CoPbase=evalin('base','CoPbase');
        IE_cop = evalin('base','IE_cop');
        weight=evalin('base','weight');
        tweight=50*0.01*weight;
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
        
        
    end
    if(get(handles.iestanding,'Value')==1)
        act = 0;
        assignin('base','TMVC',act);
        trial = 2;
        assignin('base','trial',trial);
        if act ==0
            range=1;
        else
            range=0.025
        end
        assignin('base','range',range);
        CoP=0;
        weight=evalin('base','weight');
        assignin('base','CoP',CoP);
        DP_cop = evalin('base','DP_cop');
        CoPbase=evalin('base','CoPbase');
        IE_cop = evalin('base','IE_cop');
        weight=evalin('base','weight');
        tweight=50*0.01*weight;
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
    end
    
end
if(get(handles.walking,'Value')==1)
    if(get(handles.dpwalking,'Value')==1)
        act = 0;
        assignin('base','TMVC',act);
        trial = 1;
        assignin('base','trial',trial);
        if act ==0
            range=1;
        else
            range=0.025
        end
        assignin('base','range',range);
        CoP=0;
        weight=evalin('base','weight');
        assignin('base','CoP',CoP);
        DP_cop = evalin('base','DP_cop');
        CoPbase=evalin('base','CoPbase');
        IE_cop = evalin('base','IE_cop');
        weight=evalin('base','weight');
        tweight=50*0.01*weight;
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
        
        
    end
    if(get(handles.iewalking,'Value')==1)
        act = 0;
        assignin('base','TMVC',act);
        trial = 1;
        assignin('base','trial',trial);
        if act ==0
            range=1;
        else
            range=0.025
        end
        assignin('base','range',range);
        CoP=0;
        weight=evalin('base','weight');
        assignin('base','CoP',CoP);
        DP_cop = evalin('base','DP_cop');
        CoPbase=evalin('base','CoPbase');
        IE_cop = evalin('base','IE_cop');
        weight=evalin('base','weight');
        tweight=50*0.01*weight;
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
    end
    
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


% --- Executes on button press in cop.
function cop_Callback(hObject, eventdata, handles)
% hObject    handle to cop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cop


% --- Executes on button press in standing.
function standing_Callback(hObject, eventdata, handles)
% hObject    handle to standing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of standing


% --- Executes on button press in walking.
function walking_Callback(hObject, eventdata, handles)
% hObject    handle to walking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of walking


% --- Executes when selected object is changed in uibuttongroup5.
function uibuttongroup5_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup5
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in weight.
function weight_Callback(hObject, eventdata, handles)
% hObject    handle to weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
weightcalc;
assignin('base','tweight',tweight);
assignin('base','weight',weight);
% --- Executes on button press in pbase.
function pbase_Callback(hObject, eventdata, handles)
% hObject    handle to pbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
right_loadcell_bias = pbase('PBASE.DAT');
assignin('base','f1_off',right_loadcell_bias.f1_off);
assignin('base','f2_off',right_loadcell_bias.f2_off);
assignin('base','f3_off',right_loadcell_bias.f3_off);
assignin('base','f4_off',right_loadcell_bias.f4_off);
assignin('base','f5_off',right_loadcell_bias.f5_off);
assignin('base','f6_off',right_loadcell_bias.f6_off);
assignin('base','f7_off',right_loadcell_bias.f7_off);
assignin('base','f8_off',right_loadcell_bias.f8_off);

left_loadcell_bias = pbase('PBASE_L.DAT');
assignin('base','f1_off_L',left_loadcell_bias.f1_off);
assignin('base','f2_off_L',left_loadcell_bias.f2_off);
assignin('base','f3_off_L',left_loadcell_bias.f3_off);
assignin('base','f4_off_L',left_loadcell_bias.f4_off);
assignin('base','f5_off_L',left_loadcell_bias.f5_off);
assignin('base','f6_off_L',left_loadcell_bias.f6_off);
assignin('base','f7_off_L',left_loadcell_bias.f7_off);
assignin('base','f8_off_L',left_loadcell_bias.f8_off);



function ttime_Callback(hObject, eventdata, handles)
% hObject    handle to ttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ttime as text
%        str2double(get(hObject,'String')) returns contents of ttime as a double


% --- Executes during object creation, after setting all properties.
function ttime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in visual.
function visual_Callback(hObject, eventdata, handles)
% hObject    handle to visual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of visual
