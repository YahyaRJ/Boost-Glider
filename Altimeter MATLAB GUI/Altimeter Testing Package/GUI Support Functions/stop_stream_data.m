function stop_stream_data(hObject,handles)
%STOP_STREAM_DATA Stops data stream from ASEN 2004 Altimeter Rev. 9.

read_error = [];%Create variable for later use
global isstream;%declare global variable.
if ~isempty(instrfind)%close any open serial objects
    fclose(instrfind);
end
fclose('all');%close all file IDs
s = serial(handles.port,'BaudRate',handles.baudrate);%create serial object
fopen(s);%open serial object for communications
fprintf(s,'%s\n','IDLE');%send idle command to altimeter
if(~isempty(instrfind))%close any open serial objects
    fclose(instrfind);
end
fclose('all');%close all file ID's
isstream = 0;%set stream variable to 0 to tell other funciton to stop stream
guidata(hObject,handles);%update HANDLES
end