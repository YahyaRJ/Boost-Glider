clear; clc; close all;

%% ThrustTest_Single Summary
% First, some of the summary from Thrust:
% This funciton will take in the file location of the two test setups and
% using the file names in those directories, will pull out all of the
% avialable tests, process them, and fit the data into a standard
% formatting for output
%
% ThrustTest_Single is meant to be a testing function where you can test
% the conditioning steps that your group develops. This version of the code
% loads in only one data set at a time so that the workspace is less
% cluttered. It is suggested that students make plots of their thrust data
% often throughout their development in this code section to visually check
% that their conditioning is working as desired. Once conditioning is
% working on one set of data, students are engcouraged to try other single
% data sets. Once groups are satisfied, the full funciton has the same form
% as this, so the conditioning code can simply be dragged and dropped into
% the full "Thrust.m" funciton.
% 
% Note that there is one major step missing from this testing funciton.
% That is averaging over multiple tests from the same setup. This will need
% to be added in the full "Thrust.m" function


%% Outputs:
% ThrustCurves:
%   A table containing 0.5 seconds of thrust data for each of the cases
%   available, this data will have formatting such that there are 501
%   evenly spaced thrust data points (rows) for each test (columns). The
%   ordering of the columns will go from max to min water volume in the 2L
%   bottle and then max to min in the 1.25L bottle
%
% Time:
%   A 1D array corresponding to the times of the thrust data points in the
%   ThrustCurves table
%
% <User defined variables for statistics>

% Set known sampling frequency
f= 1652; % [Hz]

%% Preallocate variables of interest
Time = 0:0.001:0.5; % just go ahead and define this, note that it will be 501 long
ThrustCurves = zeros(length(Time),1);

%% List what configuration is being used
bottleSize = 2000; % [ml]
waterSize = 900;% [ml]
% Be sure that the above matches up with the file name specified below. In
% the full "Thrust" script, all of this data will be loaded automatically
% for you
testName = 'Variable Water Volume/Test_T01_W0500_B2000'; % Should be a string of the path to the data (including the data file name)

% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////

%% Load data
% This should not have to be modified
fileName = testName; % again weird indexing is due to string arrays, we have to ask for all the characters in a row
data = readmatrix(fileName); % load the data
data = data(:,3)*4.448; % take only the third column and converting from lbf to N


%% Data Conditioning

cleaned_data = zeros(f/2,1);
step = 0:length(cleaned_data);
cleaned_time = step*1/f;
gravity = 9.81;
water_weight = waterSize/1000*gravity;

data(data(:)<4) = []; %Cleaning beginning data (4 Newtons)
for i = 0:f/2 %Using data from start to end (.5 seconds)
    cleaned_data(i+1) = data(i+1+mod(f,20));
end
cleaned_data = interp1(cleaned_time,cleaned_data,Time);
%Created indexs for max thrust and estimation of where water is expelled
[m,i] = max(cleaned_data);
[mi,ind] = min(cleaned_data(Time>.2));
ind = 200+ind; %Index correction because of the comparison used above

Actual_Thrust_Time = ind/1000

Peak_Thrust = m

%Linear equation for the change in weight over time
offset_eq = @(Time) (water_weight)/(Time(ind)*1000-...
    Time(i)*1000)*Time + water_weight;
offset_data = offset_eq(Time);
%Correcting the data with offset consideration
cleaned_data = cleaned_data-offset_data; 


standard_dev = std(cleaned_data);

Peak_Thrust_std = Peak_Thrust/standard_dev

figure(1);
plot(Time,cleaned_data);

%% Data Fitting

%Code that creates # of lines of best fits according to given stepsize
%{ 
stepsize = 20;
g = 0;
k = 0;

for i = 0:stepsize-1
    n = .5/stepsize;
    if i == stepsize -1
        k = 1;
    end
    newT = Time(Time >= i*n & Time <= (i+1)*n + k*.005);
    newD = cleaned_data(Time >= i*n & Time <= (i+1)*n+ k*.005);
    p = polyfit(newT,newD,3);
    fit = polyval(p,newT);
    figure(1);
    hold on;
    plot(newT,fit);
    for j = 1:length(fit)
        tot(j+i*length(fit))= fit(j);
    end
end
%}


%Splits the data into two portions.
%First section records max thrust and the second section creates a fit for
%the rest of the data. The two fits are then added into one congruent line
%of best fit

T1 = Time(Time <= Time(i));
T2 = Time(Time > Time(i));
D1 = cleaned_data((Time <= Time(i)));
D2 = cleaned_data((Time > Time(i)));
P1 = polyfit(T1,D1,5);
P2 = polyfit(T2,D2,7);
fit1 = polyval(P1,T1);
fit2 = polyval(P2,T2);

tot_fit = zeros(1,501);
for i = 1:length(fit1)
    tot_fit(i) = fit1(i);
end
for i = 1:length(fit2)
    tot_fit(i+length(fit1)) = fit2(i);
end
figure(1);
hold on;
plot(Time,tot_fit)



%% Sample onto the standard output array format

% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
