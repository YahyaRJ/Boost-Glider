function stream_data(hObject,handles)
% STREAM_DATA(hObject,handles) ============================================
% This function opens up the com port specified in the AltimeterGUI and
% streams pressure/altitude data depending on which radio button is
% selected within the GUI.

% First, verify that the com port is functional ===========================
if(check_port(handles.port) == 0)
    % If not then throw an error and quit the stream
    set(handles.status_text,'String', 'Com port error. Try again after changing Baud Rate');
    set(handles.status_text,'BackgroundColor', [1,0,0]);
    disp('Error: Com port is not available. Try changing Baud Rate and try again.')
    handles.status = 'Com Port Error. Try again after changing Baud Rate';
    return;
end

% Initialize variables for later use ======================================
read_error = []; % Variable to hold the error structure from a try/catch
P = []; % Pressure data array
err = 0; % Error counter
cal = zeros(1,6); % Calibration coefficients
global isstream; % Global flag for whether or not we're streaming data
isstream = 1;
% Check if the user wanted altitude or pressure
if(get(handles.feet_radiobutton,'Value'))
    feet = 1; % Altitude streaming
else
    feet = 0; % Pressure (mBar) streaming
end

% Update the GUI display ==================================================
guidata(hObject,handles);
plot(handles.axes1,P);
set(handles.axes1,'XMinorGrid','on')

% Initialize the serial port for reading ==================================
if ~isempty(instrfind)
    fclose(instrfind);
end
s = serial(handles.port,'BaudRate',handles.baudrate);
fopen(s);
if(s.BytesAvailable)
    fread(s, s.BytesAvailable);
end

% Get the calibration coefficients ========================================
fprintf(s,'%s','CONST');
for i = 1:6
    temp = fread(s,2,'uint8');
    cal(i) = bi2de([de2bi(temp(2),8) de2bi(temp(1),8)]);
end
pause(.004)

% Get a raw temperature measurement =======================================
temp = fread(s,3,'uint8');
% This value will be used for all pressure calculations while streaming
traw = bi2de([de2bi(temp(3),8) de2bi(temp(2),8) de2bi(temp(1),8)]);
pause(.01)

% Loop until 'Stop Stream' button is pressed ==============================
i = 1;
tic
P = [];
Time = [];
window_width = 15; % width of paning window, seconds
while isstream && err <= 5
    try % See if we can read a byte
        temp = fread(s,1,'uint8');
        if isempty(temp)
            err = err + 1;
        end
        if(temp == 80) % If it's a pressure reading then process it
            temp = fread(s,3); % Get the 3 bytes of data
            if(numel(temp)==3) % Make sure we got 3 bytes
                % Get the raw pressure measurement
                praw = bi2de([de2bi(temp(3),8) de2bi(temp(2),8) de2bi(temp(1),8)]);
                % Calculate the true pressure value using the calibration
                % coefficients found earlier
                P = [P CalcPres(cal(1),cal(2),cal(3),cal(4),cal(5),cal(6),praw,traw)];
                Time = [Time toc];
                % See if we've run long enough to start paning the window
                if Time(end) > window_width 
                    start_ind = find(Time >= Time(end) - window_width,1,'first');
                    plot_Time = Time(start_ind:end);
                    plot_P = P(start_ind:end);
                    xmin = plot_Time(1);
                    P = plot_P;
                    Time = plot_Time;
                else % Otherwise just plot all the data we have so far
                    plot_Time = Time;
                    plot_P = P;
                    xmin = 0;
                end
                % Zero out any pressure values that would throw warnings
%                 plot_P(plot_P < 830) = 830;
%                 plot_P(plot_P > 845) = 845;
                xmax = plot_Time(end);
                if(feet) % See if we're plotting altitude
                    % We have to convert from mBar to Pa
                    h = atmospalt(100*plot_P);
                    % If this is the first loop save the initial altitude
                    % so that everything else is relative
                    if i == 1
                        init_alt = h;
                    end
                    % Subtract the initial altitude out so that relative
                    % altitude is displayed
                    H = h - init_alt;
                    % Now convert from m to ft
                    H = H*3.28084;
                    % Plot the relative altitude
                    plot(handles.axes1,plot_Time,H);
                    ylabel('Relative Altitude, ft')
                    % Determine the vertical range of the data
                    ymin = min(H) - 0.2*range(H);
                    ymax = max(H) + 0.2*range(H)+0.1;
                else % Otherwise plot the pressure in mBar
                    plot(handles.axes1,plot_Time,plot_P);
                    ylabel('Pressure, mBar')
                    % Determine the vertical range of the data
                    ymin = min(plot_P) - 0.2*range(plot_P);
                    ymax = max(plot_P) + 0.2*range(plot_P)+0.1;
                end
                % Update the axis range
                axis([xmin xmax ymin ymax]);
                xlabel('Elapsed Time, s')
                grid on
                drawnow;
                i = i+1;
                % Clear the serial buffer
                if s.BytesAvailable > 0
                    buffer = fread(s, s.BytesAvailable);
                end
            end
        end
    catch read_error;  % If we encouter an error increase a counter and
        err = err + 1; % display the info in the command window
        if isstream % Only display error if the user hasn't stopped the stream
            disp(read_error.message)
            disp(['Loop Counter: ' num2str(i)])
        end
    end
end

if err > 5
    errordlg('ERROR -- Streaming data encountered an error and had to quit. Contact Trudy Schwartz or a Lab Assistant for help')
end
    
if ~isempty(instrfind)
    fclose(instrfind);
    delete(instrfind)
end
fclose('all');

