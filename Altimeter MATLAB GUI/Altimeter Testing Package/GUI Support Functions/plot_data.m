function plot_data(handles)
%PLOT_DATA Plots data from the filename specified by user in GUI fig. Takes
%in GUI structured variable HANDLES
%
% See also ATMOSPALT.

% Clear the display axes
cla(handles.axes1)
axes(handles.axes1)

% First try to read data to make sure file is valid
try
    filename = get(handles.processed_data_file_text,'String');%Gets filename
    data = load(filename);%load data
catch err % If it wasn't valid then throw an error but don't crash everything
    str = sprintf(['ERROR -- Invalid processed data file...\n' err.identifier]);
    errordlg(str,'Error in Plotting Data');
    return
end
    
try
    H = atmospalt(100*data(:,2));%calcualtes height based on pressure data
    H = H*3.28084; % convert from meters to feet
    if(get(handles.feet_radiobutton,'Value'))%If the feet radio button selected
        plot(data(:,1),H-H(1));%plot in feet
        ylabel('Relative Altitude, ft')
    else%if mBarr radio button is selected plot in mBarr
        plot(data(:,1),data(:,2))
        ylabel('Pressure, mBar')
    end
    xlabel('Time, s')
    xlim([data(1,1) data(end,1)])
    grid on
catch err % If it wasn't valid then throw an error but don't crash everything
    if isempty(data)
        str = ['ERROR -- tried to plot from empty data file: ' filename];
    else
        str = sprintf(['ERROR -- Invalid processed data file...\n' err.identifier]);
    end
    errordlg(str,'Error in Plotting Data');
end
