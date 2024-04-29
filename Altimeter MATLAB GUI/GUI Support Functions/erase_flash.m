function updated_handles = erase_flash(handles)
%ERASE_FLASH Sends command to erase flash to ASEN 2004 Altimeter Rev. 9
%   ERASE_FLASH(HANDLES) takes in GUI structured variable 'handles' and
%   sends the erase command over USART to the altimeter.

set(handles.status_text,'String','');%resets status string
set(handles.status_text,'BackgroundColor',[.941,.941,.941]);%resets color
if(check_port(handles.port) == 0)%standard port check...
    set(handles.status_text,'String', 'Com port error. Try again after changing Baud Rate');
    set(handles.status_text,'BackgroundColor', [1,0,0]);
    disp('Error: Com port is not available. Try changing Baud Rate and try again.')
    handles.status = 'Com Port Error';
    updated_handles = handles;
    return;
end
read_error = [];%creates variable for later use
if ~isempty(instrfindall)%closes any open serial objects
    fclose(instrfindall)
    delete(instrfindall)
end
s = serial(handles.port,'BaudRate',handles.baudrate);%creates serial object
fopen(s);%opens serial object for communications

tic
fprintf(s,'%s\n','WIPE');%sends the erase command followed by newline

disp([datestr(now) ' -- Sent WIPE command'])
fprintf('%s -- Waiting for response...',datestr(now))
ii = 1;
while s.BytesAvailable == 0 && toc < 10
    pause(0.01)
    if toc > ii
        fprintf('.')
        ii = ii + 1;
    end
end
tic
if s.BytesAvailable == 0
    fprintf('\n%s -- Got no response from altimeter!\n',datestr(now))
    handles.status = 'Unresponsive Altimeter';
else
    response = [];
    fprintf('\n%s -- Got response: ',datestr(now))
    while s.BytesAvailable ~= 0 && toc < 10
        temp = fread(s,1,'uchar');
        fprintf('%s',temp)
        response = [response char(temp)];
    end
    fprintf('\n')
    if strcmp(response, 'Memory Has Been Erased!')
        handles.status = response;
    else
        handles.status = ['Unexpected Response: ' response];
    end
end

if(~isempty(instrfindall))%closes all open serial objects
    fclose(instrfindall);
    delete(instrfindall);
end
fclose('all');%closes all file IDs
updated_handles = handles;
end