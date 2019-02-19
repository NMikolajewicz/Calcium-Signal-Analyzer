
function varargout = analyzer_UI(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyzer_UI_OpeningFcn, ...
                   'gui_OutputFcn',  @analyzer_UI_OutputFcn, ...
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

function analyzer_UI_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);

function varargout = analyzer_UI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)
file_name = get(hObject,'String');
handles.properties.file = file_name{1};
guidata(hObject,handles)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
file_name = uigetfile('*.xls*','Select excel file for preparation');
try; handles.properties; catch; handles.properties = []; end;
try; file_name = erase(file_name, '.xlsx'); catch; end;
try; file_name = erase(file_name, '.xls'); catch; end;
set(handles.edit1, 'String', file_name);
handles.properties.file = file_name;
guidata(hObject,handles);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
handles.properties.plotFigure = get(hObject, 'Value');
guidata(hObject,handles)

% --- Executes on button press in checkbox1.
function checkbox4_Callback(hObject, eventdata, handles)
handles.properties.saveFigure = get(hObject, 'Value');
guidata(hObject,handles)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global fig_save plt
try;
    whos
    assignin('base', 'hObject', hObject);
    assignin('base', 'eventdata', eventdata);
    assignin('base', 'handles', handles);
    analyzer_GUI(handles.properties);
catch e;
    disp(e)
    msgbox({'Error while running Calcium Analyzer',...
        '                                                              ',...
        'Ensure inputs are correct, refer to user guide'},'Error', 'Error');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
close all; clear all; 
msgbox({'Calcium Analyzer has been closed'}, 'Cancelled');
error('Calcium Analyzer has been closed');
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
