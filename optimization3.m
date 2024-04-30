%% Clean Workspace and Housekeeping

clear
clearvars
clc
close all

% removes warnings for table variable names for a cleaner output
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

%% Import and Read Aircraft Design File
Design_Input = readtable("Design Input File5.xlsx",'Sheet','Input','ReadRowNames',true); %Read in Aircraft Geometry File

Count = height(Design_Input); %Number of different aircraft configurations in geometry file

% Import Airfoil Data File
Airfoil = readtable("Design Input File5.xlsx",'Sheet','Airfoil_Data'); %Read in Airfoil Data

% Import Benchmark Aircraft Truth Data
Benchmark = readtable("Design Input File5.xlsx",'Sheet','Benchmark_Truth'); %Read in Benchmark "Truth" Data for model validation only

% Import Material Properties Data
Material_Data = readtable("Design Input File5.xlsx",'Sheet','Materials'); %Read in prototyp material densities for weight model

%% Optimization

%{

to meet the goal of minimum of 100m glide range, 

assume wingspan from 0-1 m
VH from 0.4-0.7
Vv from 0.04-0.07

%}

%% Caluations - Conditions and Sizing
% US Standard Atmophere - uses provided MATLAB File Exchange function
[rho,a,T,P,nu,z]= atmos(Design_Input.altitude_o(:,:)); 

ATMOS = table(rho,a,T,P,nu,z); % Reorganize atmopheric conditions into a table for ease of passing into functions
clearvars rho a T P nu z % Clear original variables now that they are in a table

% Call Wing Geometry Calcuation Function
WingGeo_Data = WingGeo(Design_Input,Count); %Calculate specific wing geometry from wing configuration parameters

%% Calculations - Lift and Drag
% Call Wing Lift & Drag Model Function
[WingLiftModel,AoA,AoA_Count,AirfoilLiftCurve,WingLiftCurve,WingDragCurve] =...
    WingLiftDrag(Design_Input,Airfoil,Count); 

% Call Parasite Drag Buildup Model Function
[Parasite_Drag_Data,FF_Table] = ...
    ParasiteDrag(Design_Input,Airfoil,WingGeo_Data,ATMOS,Count);

% Call Induced Drag Model Function
InducedDrag_Data = ...
    InducedDrag(Design_Input,WingLiftModel,WingLiftCurve,WingDragCurve,WingGeo_Data,Count,Airfoil,Parasite_Drag_Data);

% Call Complete Drag Polar Function
[DragPolar_mod1,DragPolar_mod2,DragPolar_mod3] = ...
    DragPolar(Parasite_Drag_Data,InducedDrag_Data,Design_Input,AoA_Count,WingLiftCurve,Count);

% Call L/D Analysis Function
[LD_mod1,LD_mod2,LD_mod3,LD_benchmark] = ...
    LD(Benchmark,DragPolar_mod1,DragPolar_mod2,DragPolar_mod3,WingLiftCurve,AoA_Count,Count);

% Call Weight Model
[Weight_Data,CG_Data] = ...
    Weight(Design_Input,Count,WingGeo_Data,Airfoil,Material_Data);

%% Calculations - Dynamic Models
% Call Thrust Model
[ThrustCurves, Time] = Thrust();

% Call Boost-Ascent Flight Dynamics Model
[apogee, hApogee, stateStruct] = BoostAscent(Design_Input, ATMOS, Parasite_Drag_Data, Weight_Data, ThrustCurves, Time, Count);

apogee2 = ones(Count,1) .* 16;

% Call Glide Flight Dynamics Model
[GlideRange] = GlideDescent(LD_mod1, apogee, Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count);
[GlideRange2] = GlideDescent(LD_mod1, apogee2, Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count);

%% Rough 2D Glide Range Estimation
figure()
handles = [];
labels = {};
for n=1:Count
    Range = GlideRange{n, 1}; % Range
    Range2= GlideRange2{n,1};
    height = apogee(n); %height
    height2 = apogee2(n);
    [dumby_v, idistance] = min(stateStruct.(['Config_',num2str(n)]).data(:, 6));
    rocket_distance = stateStruct.(['Config_',num2str(n)]).data(idistance,4);

    x1 = [0 Range];
    y1 = [height 0];
    pl1 = line(x1,y1);
    pl1.Color = getColor(n,Count);
    pl1.LineWidth = 1.5;

    hold on;
    x2 = [0 Range2];
    y2 = [height2 0];
    pl2 = line(x2, y2);
    pl2.Color = getColor(n,Count);
    pl2.LineWidth = 1.5;

    handles = [handles,pl1];
    labels = [labels, sprintf('Config %d',n)];
end
hold off;
xlabel('Distance (m)');
xline(100,"LineStyle","--","Label","Design Requirement");
ylabel('Height (m)');
title('Rough Glide Range Estimation');
grid on;
%xlim([-10 120]);
axis equal;
legend(handles,labels,"Location","bestoutside");




%% Finenss Ratio Optimization
figure();
fuslength = 0.29:0.1:1.19;
yyaxis left;
plot(fuslength,apogee(2:11),"LineWidth",2);
ylabel('Rocket Apogee(m)');
yyaxis right;
plot(fuslength,GlideRange{2:11, 1},"LineWidth",2);
ylabel('Glide Gange(m)');
xlabel('Fuselage Length(m)');
title('Fuselage Length Optimization');

%% Fuselage Diameter Optimization

figure();
fusdiam = 0.09:0.004:0.126;
yyaxis left;
plot(fusdiam,apogee(12:21),"LineWidth",2);
ylabel('Rocket Apogee(m)');
yyaxis right;
plot(fusdiam,GlideRange{12:21, 1},"LineWidth",2);
ylabel('Glide Gange(m)');
xlabel('Fuselage Diameter(m)');
title('Fuselage Diameter Optimization');


%% Functions

function colorscheme = getColor(index,totalCount)
    cmap = colormap('jet');
    colorIndex = round((index - 1)*(size(cmap,1)-1) / (totalCount -1)) +1;
    colorscheme = cmap(colorIndex,:);
end