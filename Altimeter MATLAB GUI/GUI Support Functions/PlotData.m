%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Function for Altimeter rev 8
% Input of filename (including .txt)
% Outputs array of broken indicies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P_Time,P] = PlotData(filename,handles)
Coefficients = [];
T_Raw = [];
P_Raw = [];
T_Time = [];
P_Time = [];
P = [];
T = [];
file = fopen(filename);

[data,count] = fread(file);
for i = 1:2:12
    temp = bi2de([de2bi(data(i+1),8) de2bi(data(i),8)]);
    Coefficients = [Coefficients temp];
end

disp(Coefficients)

while (i < count+1)
    if(data(i) == 84)           % Temperature
        temp = bi2de([de2bi(data(i+3),8) de2bi(data(i+2),8) de2bi(data(i+1),8)]);
        T_Time = [T_Time temp];
        temp = bi2de([de2bi(data(i+6),8) de2bi(data(i+5),8) de2bi(data(i+4),8)]);
        T_Raw = [T_Raw temp];
        i = i + 8;
    elseif(data(i) == 80)       % Pressure
        temp = bi2de([de2bi(data(i+3),8) de2bi(data(i+2),8) de2bi(data(i+1),8)]);
        P_Time = [P_Time temp];
        temp = bi2de([de2bi(data(i+6),8) de2bi(data(i+5),8) de2bi(data(i+4),8)]);
        P_Raw = [P_Raw temp];   %#ok<*AGROW>
        i = i + 8;
    else
        i = i + 1;
    end
end
length(P_Raw)
length(T_Raw)
for i = 1:length(P_Raw)
    temp = (T_Raw(i) + T_Raw(i+1)) / 2;
    [P(i),T(i)] = CalcPres(Coefficients(1),Coefficients(2),Coefficients(3),Coefficients(4),Coefficients(5),Coefficients(6),P_Raw(i),temp);
end

P_Time = P_Time / 1000;
if(get(handles.feet_radiobutton,'Value'))
    H = atmospalt(P*100);
    plot(P_Time,H)
else
    plot(P_Time,P);
end
fclose(file);

end