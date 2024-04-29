function process_data(handles)
%PROCESS_DATA Process raw data from ASEN 2004 Altimeter Rev. 9.
%   PROCESS_DATA(HANDLES) takes in GUI structured variable HANDLES and
%   loads file speciefied by user in GUI fig and process data to get time
%   stamps, temperature and pressure data.  Saves data to file also
%   specified by user (prompted while running).

% First try to read data to make sure file is valid
try
    filename = get(handles.raw_data_file_text,'String');%Get raw data filename
    file = fopen(filename);%open file
    [data,count] = fread(file);%read data and store in data
catch err % If it wasn't valid then throw an error but don't crash everything
    str = sprintf(['Invalid raw data file...\n' err.identifier]);
    errordlg(str,'Error in Processing Data');
    return
end

% success_flag = 0;
% while ~success_flag
%     
%     
%     [handles.process_data_filename, handles.process_data_pathname]...
%         =uiputfile('*.*');%Prompt user to create a file for processed data
%     % Check for valid input (invalid if user pressed cancel)
%     if ~ischar(handles.process_data_filename) || ~ischar(handles.process_data_pathname)
%         waitfor(errordlg('You MUST enter a valid filename to continue...'))
%     else
%         success_flag = 1;
%     end
% end

[path, name, ~] = fileparts(filename);
handles.process_data_pathname = path;
if get(handles.remove_outliers_on_pushbutton,'Value') == 1 % Removing Outliers
    handles.process_data_filename = ['\' name '_processed_smoothed'];
else
    handles.process_data_filename = ['\' name '_processed_notSmoothed'];
end



set(handles.processed_data_file_text,'String',...
    [handles.process_data_pathname handles.process_data_filename]);
%^Loads file name and path into plot data string for convenience^

wb = waitbar(0,'Processing Data... Please Wait');%Create waitbar
Coefficients = zeros(6,1);
counter = 1;
for i = 1:2:12%for the first 12 bytes
    coeff_val = bi2de([de2bi(data(i+1),8) de2bi(data(i),8)]);%convert pairs of
    %bytes to binary, combine high and low bytes and then convert back to
    %decimal
    Coefficients(counter) = coeff_val;%Each byte pair is a calibration
    %coefficient
    counter = counter + 1;
end
start_of_fname = find(filename == '\',1,'last') + 1;
name_str = ['==================== ' filename(start_of_fname:end) ' ===================='];
spacer_bar = '';
for ii = 1:length(name_str)
    spacer_bar = [spacer_bar '='];
end
disp(spacer_bar)
disp(name_str)
disp('Parsed Altimeter Calibration Coefficients:')
disp(Coefficients)%Display the coefficients, mostly for debugging...

T_Time = zeros(1,round(count/6));%Temperature time stamp
T_Raw = T_Time;%Raw temperature
P_Time = zeros(1,round(count/6));%Pressure time stamp
P_Raw = P_Time;%Raw pressure
T_counter = 1;
P_counter = 1;
end_of_file_counter = 0;
while (i < count+1)%While we haven't processed more than we read from file:
    if(data(i) == 84)           %if data(i) is an ASCII 'T' -> Temperature
        %Time stamp
        T_timestamp = bi2de([de2bi(data(i+3),8) de2bi(data(i+2),8) de2bi(data(i+1),8)]);
        T_Time(T_counter) = T_timestamp;
        %Raw Temperature
        T_val_raw = bi2de([de2bi(data(i+6),8) de2bi(data(i+5),8) de2bi(data(i+4),8)]);
        T_Raw(T_counter) = T_val_raw;
        i = i + 8;
        T_counter = T_counter + 1;
    elseif(data(i) == 80)       %if data(i) is an ASCII 'P' -> Pressure
        %Time stamp
        P_timestamp = bi2de([de2bi(data(i+3),8) de2bi(data(i+2),8) de2bi(data(i+1),8)]);
        P_Time(P_counter) = P_timestamp;
        %Raw Pressure
        P_val_raw = bi2de([de2bi(data(i+6),8) de2bi(data(i+5),8) de2bi(data(i+4),8)]);
        P_Raw(P_counter) = P_val_raw;
        i = i + 8;
        P_counter = P_counter + 1;
    else
        if data(i) == 255 & data(i-6:i-1) == 255
            % If we get 7 "empty" characters in a row, we have reached the
            % end of the file
            i = count;
        end
        i = i + 1;
    end
    try
        waitbar(i/count,wb,['Processing Data... ' num2str(i/count*100,2) '%'])%update waitbar
    catch
        wb = waitbar(i/count,['Processing Data... ' num2str(i/count*100,2) '%']);%Create waitbar
    end
end
fclose(file);%close raw data file
close(wb)%close waitbar

% Trim down the arrays to the proper size
T_Time = T_Time(1:T_counter-1);
T_Raw = T_Raw(1:T_counter-1);
P_Time = P_Time(1:P_counter-1);
P_Raw = P_Raw(1:P_counter-1);

if isempty(T_Raw) || isempty(P_Raw)
    disp('Total time of data set: N/A')
    disp(spacer_bar)
    disp(spacer_bar)
    errordlg(['File: ' filename(start_of_fname:end) ' contained no valid data...'])
    return
end

min_length = min([length(T_Raw) length(P_Raw)]);
P = zeros(1,min_length-1);%Actual Pressure
T = P;                       %Actual Temperature
for i = 1:min_length-1%calcualte actual pressure and temperature for every
    %pressure time stamp.
    T_ave = (T_Raw(i) + T_Raw(i+1)) / 2;
    [P(i),T(i)] = CalcPres(Coefficients(1),Coefficients(2),Coefficients(3),Coefficients(4),Coefficients(5),Coefficients(6),P_Raw(i),T_ave);
end

P_Time = P_Time / 1000;%convert pressure time to seconds from ms.
fprintf('Total time of data set: %.1f seconds\n',range(P_Time))
disp(spacer_bar)
disp(spacer_bar)

% If user wants to remove outliers
if get(handles.remove_outliers_on_pushbutton,'Value') == 1
    P = smoothdata(P, 'movmedian', 5);
end

handles.pressure_data = P;%Create variables in GUI structured variable
handles.temperature_data = T;
handles.time_data = P_Time;
VecLength = min([length(P) length(T) length(P_Time)]);
P_Time = P_Time(1:VecLength);
P = P(1:VecLength);
T = T(1:VecLength);
write_data = [P_Time;P;T];%Create variable to write data
%open new data file
file = fopen([handles.process_data_pathname handles.process_data_filename],'w+');
%write data to file
fprintf(file,'%% %s\n%% Time, s\t Pressure, mBar\t Temperature, degC\n',datestr(now));
count = fprintf(file,'%.3f %.4f %.4f\n',write_data);
fclose(file);%close file
% plot(handles.axes1,P_Time, P);%plot pressure data in mBar
plot_data(handles);
% figure
% plot(P_Time,T)
end