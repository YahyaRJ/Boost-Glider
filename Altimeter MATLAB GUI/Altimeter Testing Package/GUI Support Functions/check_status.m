function status = check_status(handles)
%CHECK_STATUS checkst the status of the ASEN 2004 Altimeter Rev. 9.
%   STATUS = CHECK_STATUS(HANDLES) takes the GUI structured variable
%   'HANDLES' and communicates with the altimeter to see what its current
%   status is.  Outputs STATUS and updates handles.status_text to display
%   on GUI fig.
%   STATUS is a string:
%       'launch'-Flight has been written to flash
%       'erased'-No flight has been written to flash

%Performs port check and displays error if necessary.
if(check_port(handles.port) == 0)
    set(handles.status_text,'String', 'Com port error. Try again after changing Baud Rate');%Displayes string
    set(handles.status_text,'BackgroundColor', [1,0,0]);%changes color to red
    disp('Error: Com port is not available. Try changing Baud Rate and try again.')
    %displays error message to command terminal
    status = 'Com Port Error';%changes status
    return;
end


read_error = [];%creates variable for later use
if ~isempty(instrfind)%if serial object is open,
    fclose(instrfind) %close it...
end
s = serial(handles.port,'BaudRate',handles.baudrate);%creates new serial object
fopen(s);%opens serial object for communications
fprintf(s,'%s\n','STAT');%Sends the string 'STAT' followed by a new line
                         %over USART to cause the altimeter to respond with
                         %status
                         
fprintf('%s -- Requesting Altimeter Status...\n',datestr(now))
try
    alt_status = fread(s,1,'uchar');%read 1 unsigned character.
catch read_error;%If there is an error, store it in 'read_error'
    set(handles.status_text,'String', 'Serial read error.');%Display error
    set(handles.status_text,'BackgroundColor', [1,0,0]);%Change color to red
    status = 'error';%change status
    fclose(s);%close serial object
    clear s;%clear variable for later use
    return;
end

if(alt_status == 'L')% If launch has occurred and memory contains data
    set(handles.status_text, 'String', 'Flight written to flash.');
    set(handles.status_text, 'BackgroundColor', [0,1,0]);
    status = 'launched';%change status
    disp([datestr(now) ' -- Current Altimeter Status: LAUNCHED (' sprintf('%c',alt_status) ')'])
elseif(alt_status == 'E')%If memory is erased and contains no relavant data
    set(handles.status_text, 'String', 'Flash erased and ready for flight.');
    set(handles.status_text, 'BackgroundColor', [0.6,0.6,1]);
    status = 'erased';%change status
    disp([datestr(now) ' -- Current Altimeter Status: ERASED (' sprintf('%c',alt_status) ')'])
else%If the returned status is neither an L nor E
    if isempty(alt_status)
        status_str = ['Unknown response from altimeter: ' sprintf('%c',alt_status) ' (check board connections)'];
    else
        status_str = ['Unknown response from altimeter: ' sprintf('%c',alt_status)];
    end
    set(handles.status_text, 'String', status_str);
    set(handles.status_text, 'BackgroundColor', [1,0.5,0.3]);
    status = 'unknown';%change status
    disp([datestr(now) ' -- Current Altimeter Status: UNKNOWN (' sprintf('%c',alt_status) ')'])
end
fclose('all');%Close all open file ID's (including serial objects)
if(~isempty(instrfind))%If any are still open, close them...
    fclose(instrfind);
    delete(instrfind);
end
end