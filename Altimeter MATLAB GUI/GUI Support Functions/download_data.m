function download_data(handles)
%DOWNLOAD_DATA Downloads data from ASEN 2004 Altimeter Rev. 9 and saves raw
%data file.
%   DOWNLOAD_DATA(HANDLES) uses GUI structured variable 'handles' to
%   communicate with the altimeter and retrieve data from altimeter's flash
%   memory chip.  Data is written to a file.  User is prompted for file
%   name.

%same old com port check...see check_port.m
if(check_port(handles.port) == 0)
    set(handles.status_text,'String', 'Com port error. Try changing Baud Rate and try again');
    set(handles.status_text,'BackgroundColor', [1,0,0]);
    disp('Error: Com port is not available. Change Baud Rate and try again.')
    handles.status = 'Com Port Error';
    return;
end

read_error = [];%creates variable for later use
data = [];%creates variable for later use
dataLength = length(data);
ii = 1;%creates variable for later use
wb = waitbar(0,'Downloading... Do not press any buttons...');%creates a
                                                             %waitbar
err = 0;%creates variable for later use
if ~isempty(instrfind)%Closes any open serial objects
    fclose(instrfind);
end
s = serial(handles.port,'BaudRate',handles.baudrate,'InputBufferSize',512000);%creates serial object
fopen(s);%opens port for communication
fprintf('%s -- Sent DOWNLOAD command to altimeter...\n',datestr(now))
fprintf(s,'%s\n','DOWN');%sends the download command over USART to altimeter
full_flight_flag = 0; % Flag stating whether or not the full file was successfully read from flash
tic
data_timeout = toc;
while(err < 2) && (data_timeout < 5)
% if there have been less than 2 errors (typically no errors)
% and if we get data more frequently than every 5 seconds
    
    if s.BytesAvailable > 0
        tic
    try
        data = [data; fread(s,s.BytesAvailable)];%read the buffer and store it to 'data'
    catch read_error;%If there was an error, count it and store it to var.
        % This only gets incremented if we run out of data to
        % read and we never the 'EEEEE' string
        err = err + 1;
    end
    dataLength = length(data);
    % Check to see if we get the 5 'E' characters signifying that we are at
    % the end of a full data set
    if (dataLength >=5 & data(dataLength) == 69 & data(dataLength-4:dataLength-1) == 69)
        full_flight_flag = 1;
        fprintf('%s -- Successfully downloaded full flight data!\n',datestr(now))
        break
    end
    % Check to see if we are now reading "empty" memory where the bytes are
    % just set to high
    if (dataLength >= 7 & data(dataLength) == 255 & data(dataLength-6:dataLength-1) == 255)
        full_flight_flag = 0;
        fprintf('%s -- WARNING: flight data ended before full flight cycle...\n',datestr(now))
        break
    end
    end
    data_timeout = toc;
    ii = ii+1;%increment i
    waitbar(mod(ii,5000)/5000,wb);%increase waitbar progress (each file should be
                        %about <50000 bytes so waitbar is fairly accurate
end
if(~isempty(instrfind))%closes any open serial objects
    fclose(instrfind);
    delete(instrfind);
end
close(wb);%closes waitbar

% If the full flight was not written to flash, the user must power cycle
% the device in order to get it to stop dumping its memory to the serial
% port
if ~full_flight_flag && ~isempty(data)
    alt_reset_flag = 0;
    reset_try_counter = 0;
    while ~alt_reset_flag
        % Tell the user to power cycle
        waitfor(warndlg(['WARNING: It was detected that the flight data' ...
            ' written to memory ended before a full flight cycle completed.'...
            ' Power may have been lost during flight so some data may be missing.'...
            ' The altimeter must now be power-cycled... (turn it off then on again)'],...
            'Power Cycle Altimeter'));
        
        % Check to see if the altimeter is still dumping memory
        status = check_status(handles);
        % If it's an unknown response it's probably still dumping
        if strcmp(status,'unknown')
            alt_reset_flag = 0;
            reset_try_counter = reset_try_counter + 1;
            if reset_try_counter >= 3
                waitfor(errordlg('The altimeter has not been successfully power-cycled after 3 tries and is still dumping its memory. Please see Trudy or Bobby!'));
                break
            end
        else
            alt_reset_flag = 1;
        end
    end
end

% Only perform the following if we read any data successfully =============
if ~isempty(data)
    success_flag = 0;
    while ~success_flag
        [handles.retrieve_data_filename, handles.retrieve_data_pathname]...
            = uiputfile('*.*','Create Raw Data File');%prompts user to create file
        
        % Check for valid input (invalid if user pressed cancel)
        if ~ischar(handles.retrieve_data_filename) || ~ischar(handles.retrieve_data_pathname)
            button = questdlg('You will LOSE THE DATA unless you enter a valid filename! Continue?','Lose the data?','LOSE THE DATA','Cancel','LOSE THE DATA');
            if strcmp(button,'LOSE THE DATA')
                break
            end
        else
            success_flag = 1;
        end
    end
    
    if success_flag
        file = fopen([handles.retrieve_data_pathname handles.retrieve_data_filename],'w');
        %^opens file for writing^
        count = fwrite(file,data(1:end));%writes data to file
        fprintf('%s -- Successfully saved raw data file: %s\n',datestr(now),handles.retrieve_data_filename)
        set(handles.raw_data_file_text,'String',[handles.retrieve_data_pathname handles.retrieve_data_filename]);
        %^Loads file path and name into the process data file field on GUI for
        %convenience
    end
else
    waitfor(errordlg('Failed to successfully read altimeter data...'))
    fprintf('%s -- ERROR: failed to successfully read altimeter data...\n',datestr(now))
end
fclose('all');%closes all objects just in case...
