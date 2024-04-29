function launch(handles)
%LAUNCH Sends launch command to ASEN 2004 Altimeter Rev. 9.
%   LAUNCH(HANDLES) takes GUI structured variable 'HANDLES' and sends the
%   launch command over USART to altimeter.

set(handles.status_text,'String','');%Resets status string
set(handles.status_text,'BackgroundColor',[.941,.941,.941]);%resets color
if(check_port(handles.port) == 0)%standard port check...
    set(handles.status_text,'String', 'Com port error.');
    set(handles.status_text,'BackgroundColor', [1,0,0]);
    disp('Error: Com port is not available.  Restart GUI and try again.')
    handles.status = 'Com Port Error';
    return;
end
read_error = [];%create var for later use
i = 0;
wb = waitbar(i,'Please wait for launch sequence to complete...');%create 
                                                                 %Waitbar
if ~isempty(instrfind)%close any serial objects
    fclose(instrfind)
end

response = [];
s = serial(handles.port,'BaudRate',handles.baudrate);%Create serial object
fopen(s);%open serial object for com
fprintf(s,'%s','LAUNC');%Send launch command to altimeter of USART
for i = 1:400%wait for launch sequence to finish and update waitbar
    pause(.125/2);
    waitbar(i/400,wb);
    if s.BytesAvailable ~= 0
        response_char = fread(s,1);
        response = [response char(response_char)];
        if length(response) >= length('Memory not erased, cannot Launch.')
        if strcmp(response(end-32:end),'Memory not erased, cannot Launch.')
            waitfor(errordlg('ERROR: The altimeter was unable to launch because flash memory was not cleared!'))
            break
        end
        end
    end
end
close(wb)%close waitbar.
fclose('all');
if ~isempty(instrfind)%close any serial objects
    fclose(instrfind);
end
end