%% Clean Workspace and Housekeeping

clear
clearvars
clc
close all

% removes warnings for table variable names for a cleaner output
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

%% Import and Read Aircraft Design File
Design_Input = readtable("Design Input File7.xlsx",'Sheet','Input','ReadRowNames',true); %Read in Aircraft Geometry File

Count = height(Design_Input); %Number of different aircraft configurations in geometry file

% Import Airfoil Data File
Airfoil = readtable("Design Input File7.xlsx",'Sheet','Airfoil_Data'); %Read in Airfoil Data

% Import Benchmark Aircraft Truth Data
Benchmark = readtable("Design Input File7.xlsx",'Sheet','Benchmark_Truth'); %Read in Benchmark "Truth" Data for model validation only

% Import Material Properties Data
Material_Data = readtable("Design Input File7.xlsx",'Sheet','Materials'); %Read in prototyp material densities for weight model

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

building = ones(101, 1) * 16;

% Call Glide Flight Dynamics Model
%[GlideRange] = GlideDescent(LD_mod1, [16,16,16,16], Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count);
[GlideRange] = GlideDescent(LD_mod1, apogee, Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count);
[GlideRange2] = GlideDescent(LD_mod1, building, Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count);
%% Plotting
figure
hold on
for n = 1:Count
    plot(WingLiftCurve{1,:},DragPolar_mod1{n,:}); % brace indexing for plotting tables
end
xlabel('Coefficient of Lift (CL)');
ylabel('Coefficient of Drag (CD)');
title('Drag Polar Configuration Comparison');
legend(Design_Input.Config(:),'Location','bestoutside');
hold off

%% Boost_Ascent Flight Profile Plots
% Some setup to make plots more readable in color, look up the
% documentation for 'cmap' for other color map options
fields = fieldnames(stateStruct);
figure()

cmap = colormap(jet(Count+1));
set(0,'DefaultAxesColorOrder',cmap)
set(gca(),'ColorOrder',cmap);

for n = 1:Count
    distBoost = vecnorm([stateStruct.(fields{n}).data(:, 4), stateStruct.(fields{n}).data(:, 5)], 2, 2);
    plot(distBoost,...
            -stateStruct.(fields{n}).data(:, 6), ...
            DisplayName=Design_Input.Properties.RowNames{n}, Color=cmap(n, :))
    if n == 1
        hold on
    end
end
xlabel('Total Distance Traveled [m]');
ylabel('Height Achieved [m]');
title('Boost 2D Total Distance Traveled');
legend();
grid on
hold off

figure()
for n = 1:Count
    plot(stateStruct.(fields{n}).data(:, 4),...
            stateStruct.(fields{n}).data(:, 5), ...
            DisplayName=Design_Input.Properties.RowNames{n}, Color=cmap(n, :))
    if n == 1
        hold on
    end
end
xlabel('y [m] - Positive = East');
ylabel('x [m] - Positive = North');
title('Boost Ground Track');
legend();
grid on
hold off

figure()
for n = 1:Count    
    plot3(stateStruct.(fields{n}).data(:, 4),...
            stateStruct.(fields{n}).data(:, 5),...
            stateStruct.(fields{n}).data(:, 6), ...
            DisplayName=Design_Input.Properties.RowNames{n})
    if n == 1
        hold on
    end
end
xlabel('x [m] - Positive = North');
ylabel('y [m] - Positive = East');
zlabel('z [m] - Positive = Down');
title('Boost Trajectory Plots');
set(gca, 'ZDir','reverse')
set(gca, 'YDir','reverse')
legend();
grid on
axis equal
hold off

figure()
for n = 1:Count
    Wx = -Design_Input.V_wind(n)*cosd(Design_Input.Wind_Az(n)); 
    Wy = -Design_Input.V_wind(n)*sind(Design_Input.Wind_Az(n));
    quiver3(stateStruct.(fields{n}).data(:, 4), ... % x
        stateStruct.(fields{n}).data(:, 5), ... % y
        stateStruct.(fields{n}).data(:, 6), ... % z
        stateStruct.(fields{n}).data(:, 1)-Wx, ... % Vax
        stateStruct.(fields{n}).data(:, 2)-Wy, ... % Vay
        stateStruct.(fields{n}).data(:, 3), ... % Vz
        DisplayName=Design_Input.Properties.RowNames{n}, ...
        LineWidth=2)
    if n == 1
        hold on
    end
end
xlabel('x [m] - Positive = North');
ylabel('y [m] - Positive = East');
zlabel('z [m] - Positive = Down');
title('Boost Trajectory Plots with Heading Vetors');
set(gca, 'ZDir','reverse')
set(gca, 'YDir','reverse')
legend();
grid on
axis equal
hold off

%% Reset default color order
set(0,'DefaultAxesColorOrder','default')

%% Rough 2D Glide Range Estimation
figure()
handles = [];
labels = {};
for n=1:Count
    Range = GlideRange{n, 1}; % Range
    Range2= GlideRange2{n,1};
    height = apogee(n); %height
    height2 = building(n);
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



%% Fineness Ratio

figure();
fineness = 0.4:0.012:0.988;
yyaxis left;
plot(fineness,apogee(2:51),"LineWidth",2);
ylabel('Rocket Apogee(m)');
yyaxis right;
plot(fineness,GlideRange{2:51, 1},"LineWidth",2);
ylabel('Glide Gange(m)');
xlabel('Fineness Ratio');
title('Fineness Ratio Optimization');

%% Sweep Angle 

figure();
angle = 1:1:50;
yyaxis left;
plot(angle,apogee(52:101),"LineWidth",2);
ylabel('Rocket Apogee(m)');
yyaxis right;
plot(angle,GlideRange{52:101, 1},"LineWidth",2);
ylabel('Glide Gange(m)');
xlabel('Sweep Angle(deg)');
title('Sweep Angle Optimization');

%% Functions

function colorscheme = getColor(index,totalCount)
    cmap = colormap('jet');
    colorIndex = round((index - 1)*(size(cmap,1)-1) / (totalCount -1)) +1;
    colorscheme = cmap(colorIndex,:);
end