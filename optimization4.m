%% Clean Workspace and Housekeeping

clear
clearvars
clc
close all

% removes warnings for table variable names for a cleaner output
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

%% Import and Read Aircraft Design File
Design_Input = readtable("Design Input File6.xlsx",'Sheet','Input','ReadRowNames',true); %Read in Aircraft Geometry File

Count = height(Design_Input); %Number of different aircraft configurations in geometry file

% Import Airfoil Data File
Airfoil = readtable("Design Input File6.xlsx",'Sheet','Airfoil_Data'); %Read in Airfoil Data

% Import Benchmark Aircraft Truth Data
Benchmark = readtable("Design Input File6.xlsx",'Sheet','Benchmark_Truth'); %Read in Benchmark "Truth" Data for model validation only

% Import Material Properties Data
Material_Data = readtable("Design Input File6.xlsx",'Sheet','Materials'); %Read in prototyp material densities for weight model

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


% Call Glide Flight Dynamics Model
[GlideRange] = GlideDescent(LD_mod1, apogee, Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count);
