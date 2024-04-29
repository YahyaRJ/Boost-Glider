function varargout = AltimeterGUI(varargin)
% ALTIMETERGUI MATLAB code for AltimeterGUI.fig
%      ALTIMETERGUI, by itself, creates a new ALTIMETERGUI or raises the existing
%      singleton*.
%
%      H = ALTIMETERGUI returns the handle to a new ALTIMETERGUI or the handle to
%      the existing singleton*.
%
%      ALTIMETERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALTIMETERGUI.M with the given input arguments.
%
%      ALTIMETERGUI('Property','Value',...) creates a new ALTIMETERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AltimeterGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AltimeterGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AltimeterGUI

% Last Modified by GUIDE v2.5 24-Apr-2018 14:39:16

cd GUI' Support Functions'\ %change to the gui function directory
gui_dir = pwd;      %set a variable to the gui function directory
cd ..               %change directory back to the original directory
addpath(gui_dir);   % add the gui function directory to the matlab search path
%to here

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AltimeterGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AltimeterGUI_OutputFcn, ...
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
fclose('all');%Close all open serial connections
if(~isempty(instrfind))
    fclose(instrfind);%if there are any serial objects, close them.... just
    delete(instrfind);%to be sure...
end                   

% --- Executes just before AltimeterGUI is made visible.
function AltimeterGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AltimeterGUI (see VARARGIN)

% Choose default command line output for AltimeterGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AltimeterGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = AltimeterGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in com_port_popupmenu.
function com_port_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to com_port_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns com_port_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from com_port_popupmenu


% --- Executes during object creation, after setting all properties.
function com_port_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_port_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
port = get_highest_port(); % Gets the highest available com port
set(hObject,'String',port); % Puts that port into the drop down menu on GUI
set(hObject,'Value',1); % Puts a value cause it needs one...Doesn't matter...
handles.port = port; % Creates a variable that all the functions can see.

if strcmp(port,'ERR')
    waitfor(warndlg('WARNING: No altimeter detected on startup! You will still be able to process data.'))
end

guidata(hObject, handles);%Updates the GUI structured variable 'handles'


% --- Executes on selection change in baud_rate_popupmenu.
function baud_rate_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to baud_rate_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String')); % returns baud_rate_popupmenu contents as cell array
baudrate = str2double(contents{get(hObject,'Value')}); % returns selected item from baud_rate_popupmenu
handles.baudrate = baudrate;% saves the value
guidata(hObject,handles);%updates 'handles'


% --- Executes during object creation, after setting all properties.
function baud_rate_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baud_rate_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% NOTE: The following line was removed temporarily on 1/14/2015 by Thomas
% Green as the 19200 and 57600 baud rates are not currently supported by
% the altimeters. This line can be uncommented once that functionality is
% added.
set(hObject,'String',{'9600';'19200';'57600';'115200'});
% set(hObject,'String','9600')
handles.baudrate = 9600;
guidata(hObject,handles);



function retrieve_data_filename_Callback(hObject, eventdata, handles)
% hObject    handle to retrieve_data_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of retrieve_data_filename as text
%        str2double(get(hObject,'String')) returns contents of retrieve_data_filename as a double
handles.retrieve_data_filename = get(hObject,'String');%Get the filename of
                                                       %the raw data file

guidata(hObject, handles);%updates 'handles'

% --- Executes during object creation, after setting all properties.
function retrieve_data_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to retrieve_data_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in retrieve_data_pushbutton.
function retrieve_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to retrieve_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable all buttons while data is downloaded
set(handles.check_status_pushbutton,'Enable','off')
set(handles.plot_data_pushbutton,'Enable','off')
set(handles.choose_file_2_pushbutton,'Enable','off')
set(handles.choose_file_pushbutton,'Enable','off')
set(handles.process_data_pushbutton,'Enable','off')
set(handles.stream_data_pushbutton,'Enable','off')
set(handles.stop_stream_pushbutton,'Enable','off')
set(handles.launch_pushbutton,'Enable','off')
set(handles.erase_pushbutton,'Enable','off')
set(handles.retrieve_data_pushbutton,'Enable','off')
set(handles.com_port_popupmenu,'Enable','off')
set(handles.baud_rate_popupmenu,'Enable','off')
set(handles.remove_outliers_on_pushbutton,'Enable','off')
set(handles.radiobutton5,'Enable','off')
set(handles.feet_radiobutton,'Enable','off')
set(handles.mBar_radiobutton,'Enable','off')
guidata(hObject, handles);%updates 'handles'

download_data(handles);%Downloads raw data from altimeter

% Re-enable all buttons after data is downloaded
set(handles.check_status_pushbutton,'Enable','on')
set(handles.plot_data_pushbutton,'Enable','on')
set(handles.choose_file_2_pushbutton,'Enable','on')
set(handles.choose_file_pushbutton,'Enable','on')
set(handles.process_data_pushbutton,'Enable','on')
set(handles.stream_data_pushbutton,'Enable','on')
set(handles.stop_stream_pushbutton,'Enable','on')
set(handles.launch_pushbutton,'Enable','on')
set(handles.erase_pushbutton,'Enable','on')
set(handles.retrieve_data_pushbutton,'Enable','on')
set(handles.com_port_popupmenu,'Enable','on')
set(handles.baud_rate_popupmenu,'Enable','on')
set(handles.remove_outliers_on_pushbutton,'Enable','on')
set(handles.radiobutton5,'Enable','on')
set(handles.feet_radiobutton,'Enable','on')
set(handles.mBar_radiobutton,'Enable','on')
guidata(hObject, handles);%updates 'handles'

% --- Executes on button press in check_status_pushbutton.
function check_status_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to check_status_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.status = check_status(handles);%Checks altimeter status
guidata(hObject, handles);%updates 'handles'


% --- Executes on button press in erase_pushbutton.
function erase_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to erase_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
erase_flash(handles);%Erases flash memory on altimeter
guidata(hObject, handles);%updates 'handles'
handles.status = check_status(handles);%checks altimeter status
guidata(hObject, handles);%updates 'status'... cause you cant do that too
                          %often

% --- Executes on button press in launch_pushbutton.
function launch_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to launch_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
launch(handles);%Causes altimter to launch as if RBF pin is pulled
guidata(hObject, handles);%updates 'handles'


% --- Executes on button press in stream_data_pushbutton.
function stream_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to stream_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
% Disable buttons that shouldn't be pressed while streaming
set(handles.check_status_pushbutton,'Enable','off')
set(handles.plot_data_pushbutton,'Enable','off')
set(handles.choose_file_2_pushbutton,'Enable','off')
set(handles.choose_file_pushbutton,'Enable','off')
set(handles.process_data_pushbutton,'Enable','off')
set(handles.stream_data_pushbutton,'Enable','off')
set(handles.launch_pushbutton,'Enable','off')
set(handles.erase_pushbutton,'Enable','off')
set(handles.retrieve_data_pushbutton,'Enable','off')
set(handles.com_port_popupmenu,'Enable','off')
set(handles.baud_rate_popupmenu,'Enable','off')
set(handles.remove_outliers_on_pushbutton,'Enable','off')
set(handles.radiobutton5,'Enable','off')
set(handles.feet_radiobutton,'Enable','off')
set(handles.mBar_radiobutton,'Enable','off')

stream_data(hObject,handles);%streams live pressure data from altimeter
guidata(hObject, handles);%updates 'handles'
catch error_msg
   waitfor(errordlg(['ERROR: Altimeter has failed to communicate.'...
       ' Check connections and power.']))
   stop_stream_pushbutton_Callback(handles.stop_stream_pushbutton,[],handles)
end

% --- Executes on button press in stop_stream_pushbutton.
function stop_stream_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to stop_stream_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Re-enable buttons that were disabled while streaming
guidata(hObject, handles);%updates 'handles'
set(handles.check_status_pushbutton,'Enable','on')
set(handles.plot_data_pushbutton,'Enable','on')
set(handles.choose_file_2_pushbutton,'Enable','on')
set(handles.choose_file_pushbutton,'Enable','on')
set(handles.process_data_pushbutton,'Enable','on')
set(handles.stream_data_pushbutton,'Enable','on')
set(handles.launch_pushbutton,'Enable','on')
set(handles.erase_pushbutton,'Enable','on')
set(handles.retrieve_data_pushbutton,'Enable','on')
set(handles.com_port_popupmenu,'Enable','on')
set(handles.baud_rate_popupmenu,'Enable','on')
set(handles.remove_outliers_on_pushbutton,'Enable','on')
set(handles.radiobutton5,'Enable','on')
set(handles.feet_radiobutton,'Enable','on')
set(handles.mBar_radiobutton,'Enable','on')

stop_stream_data(hObject,handles);%stops the live stream
guidata(hObject, handles);%'udpates 'handles'


% --- Executes on button press in choose_file_pushbutton.
function choose_file_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to choose_file_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename, Pathname] = uigetfile('*.*');%Opens file dialog for user to
                                        %choose a raw data file
% Check for valid input (invalid if user pressed cancel)
if ~ischar(Filename) || ~ischar(Pathname)
    Pathname = sprintf('PLEASE SELECT A VALID FILE FOR PROCESSING\n Please press the ''Choose File'' button');
    Filename = '';
end
                                        
handles.data_file = [Pathname Filename];%creates the path and name
set(handles.raw_data_file_text,'String',handles.data_file);%displays it
guidata(hObject, handles);%updates 'handles'


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in process_data_pushbutton.
function process_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to process_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable all buttons while data is processed
set(handles.check_status_pushbutton,'Enable','off')
set(handles.plot_data_pushbutton,'Enable','off')
set(handles.choose_file_2_pushbutton,'Enable','off')
set(handles.choose_file_pushbutton,'Enable','off')
set(handles.process_data_pushbutton,'Enable','off')
set(handles.stream_data_pushbutton,'Enable','off')
set(handles.stop_stream_pushbutton,'Enable','off')
set(handles.launch_pushbutton,'Enable','off')
set(handles.erase_pushbutton,'Enable','off')
set(handles.retrieve_data_pushbutton,'Enable','off')
set(handles.com_port_popupmenu,'Enable','off')
set(handles.baud_rate_popupmenu,'Enable','off')
set(handles.remove_outliers_on_pushbutton,'Enable','off')
set(handles.radiobutton5,'Enable','off')
set(handles.feet_radiobutton,'Enable','off')
set(handles.mBar_radiobutton,'Enable','off')
guidata(hObject, handles);%updates 'handles'

process_data(handles);%Processes raw data

% Re-enable all buttons after data is processed
set(handles.check_status_pushbutton,'Enable','on')
set(handles.plot_data_pushbutton,'Enable','on')
set(handles.choose_file_2_pushbutton,'Enable','on')
set(handles.choose_file_pushbutton,'Enable','on')
set(handles.process_data_pushbutton,'Enable','on')
set(handles.stream_data_pushbutton,'Enable','on')
set(handles.stop_stream_pushbutton,'Enable','on')
set(handles.launch_pushbutton,'Enable','on')
set(handles.erase_pushbutton,'Enable','on')
set(handles.retrieve_data_pushbutton,'Enable','on')
set(handles.com_port_popupmenu,'Enable','on')
set(handles.baud_rate_popupmenu,'Enable','on')
set(handles.remove_outliers_on_pushbutton,'Enable','on')
set(handles.radiobutton5,'Enable','on')
set(handles.feet_radiobutton,'Enable','on')
set(handles.mBar_radiobutton,'Enable','on')
guidata(hObject, handles);%updates 'handles'

% --- Executes on button press in plot_data_pushbutton.
function plot_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plot_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_data(handles);%plots processed data
guidata(hObject,handles);%updates 'handles'

% --- Executes on button press in choose_file_2_pushbutton.
function choose_file_2_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to choose_file_2_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename, Pathname] = uigetfile('*.*');%opens dialog for user to choose file

% Check for valid input (invalid if user pressed cancel)
if ~ischar(Filename) || ~ischar(Pathname)
    Pathname = sprintf('PLEASE SELECT A VALID FILE FOR PLOTTING\n Please press the ''Choose File'' button');
    Filename = '';
end

handles.processed_data_file = [Pathname Filename];%creates path and filename
set(handles.processed_data_file_text,'String',handles.processed_data_file);
guidata(hObject, handles);%upadates 'handles'
