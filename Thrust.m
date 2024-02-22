function [ThrustCurves, Time] = Thrust()
clc; clear; close all;
%% Thrust Summary
% This funciton will take in the file location of the two test setups and
% using the file names in those directories, will pull out all of the
% avialable tests, cendition their data, and fit that data into a standard
% formatting for output. Note that statistics are also requested for the
% student deliverable, but how to pass those out will be left up to the
% students as they are not needed to be passed into any later functions.
% Despite this, the first two outputs of the funciton are not permitted to
% have their form modified.

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
% <User defined variable(s) for statistics>
%

%% Define data locations
% This is hard coded!!!
fileLoc_2L = 'Seperated Data/Variable Bottle Volume 60 psi/2000mL/'; % path to the data files, be sure to include a trailing slash
fileLoc_1pt25L = 'Seperated Data/Variable Bottle Volume 60 psi/1250mL/'; % path to the data files, be sure to include a trailing slash

%% Read in all of the avilable data and find what data there is
testInfo_2L = getThrustTestNames(fileLoc_2L);
configs_2L = unique(testInfo_2L.waterVol);
numConfigs_2L = length(configs_2L);


testInfo_1pt25L = getThrustTestNames(fileLoc_1pt25L);
configs_1pt25L = unique(testInfo_1pt25L.waterVol);
numConfigs_1pt25L = length(configs_1pt25L);


numConfigs = numConfigs_2L + numConfigs_1pt25L;

% Set known sampling frequency
f= 1652; % [Hz]

%% Preallocate variables of interest
Time = 0:0.001:0.5; % just go ahead and define this, note that it will be 501 long
ThrustCurves = zeros(length(Time),numConfigs);

ThrustCurvesNames = {};

%% Loop over all of the configurations
for N = 1:numConfigs % use upper case N to distiguish that it is counting something different from the aerodynamic modeling loops
    %% Dertemine what configuration to use for this iteration in the loop
    if N <=  numConfigs_2L % determine if we should be reading 2L or 1.25L data
        bottleSize = '2000'; % [ml]
        waterSize = configs_2L(N);
        testIndexes = find(testInfo_2L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_2L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    else
        bottleSize = '1250'; % [ml]
        waterSize = configs_1pt25L(N-numConfigs_2L);
        testIndexes = find(testInfo_1pt25L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_1pt25L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    end

    % /////////////////////////////////////////////////////////////////////////
    % MODIFY THIS SECTION
    % /////////////////////////////////////////////////////////////////////////
    % Notice that there is little to no guidance in place for this
    % function. This is on purpose as there are many different and equally
    % valid ways to process data (not to say that any way is valid though).
    % The lack of guidance is therefore to encourage you to think about,
    % discuss, and potentially debate as a group the best set of steps to
    % extract just the meaningful part of the thrust profile
    tot_data = zeros(1,501);
    num = 0;
    for i = 1:numTests
    %% Load data
        % The folloowing three lines will pull all of the files in each
        % test setup for you and give the array "data" which should be
        % conditioned. You should not need to modify any of this section of
        % code
        fileName = testNames(i, :); % again weird indexing is due to string arrays, we have to ask for all the characters in a row
        data = readmatrix(fileName); % load the data
        data = data(:,3)*4.448; % take only the third column and converting from lbf to N

        %% Data Conditioning
        data(isnan(data)) = [];
        
        if mean(data) > 0
            
        cleaned_data = zeros(f/2,1);
        step = 0:length(cleaned_data);
        cleaned_time = step*1/f;
        gravity = 9.81;
        water_weight = waterSize/1000*gravity;
        
        
        data(find(data)<10) = [];
        index = find(data > 25);
        
        
        %data(data(:)<3) = []; %Cleaning beginning data (4 Newtons)
        for j = 0:f/2 %Using data from start to end (.5 seconds)
            cleaned_data(j+1) = data(j+1+index(1)-50); %if something breaks try doing data(j+1+mod(f,20)) instead of data(j+1)
        end
        
        
        %{
        data = data(~isnan(data));

        stdData = std(data);
        newData = [];
        deviationDetected = false;
        for j = 1:length(data)
            if ~deviationDetected && abs(data(j) - mean(data)) > (1*stdData)
                deviationDetected = true;
            end
            if deviationDetected
                newData = [newData; data(j)];
            if length(newData) >= length(Time)
                break;
            end
            end
        end
        %}

        
        cleaned_data = interp1(cleaned_time,cleaned_data,Time);
        %Created indexs for max thrust and estimation of where water is expelled
        [m,in] = max(cleaned_data);
        [mi,ind] = min(cleaned_data(Time>.20));
        ind = 200+ind; %Index correction because of the comparison used above
        endData  = cleaned_data(Time>.45);
   
        Actual_Thrust_Time(i) = ind/1000;
        
        Peak_Thrust(i) = m;
        water_weight = mean(endData);
        %Linear equation for the change in weight over time
        offset_eq = @(Time) (water_weight)/(Time(ind)*1000-...
            Time(in)*1000)*Time + water_weight;
        offset_data = offset_eq(Time);
        offset_data(Time < Time(in)) = 0;
        %Correcting the data with offset consideration
        cleaned_data = cleaned_data-offset_data; 
        
        
        standard_dev = std(cleaned_data);
        
        %{
        figure(N)
        hold on;
        plot(Time,cleaned_data)
        %}
        

        
        
    %% Averaging

    dataN(i,:) = cleaned_data;
    uncertainty(i) = std(cleaned_data);

    num = num +1;
        end
    end
    for i = 1:num
        tot_data = tot_data + dataN(i,:);
    end
    
    tot_data = tot_data./num;


    std_thrustTime(N) = std(Actual_Thrust_Time)
    avg_thrustTime(N) = mean(Actual_Thrust_Time)
    
    avg_maxThrust(N) = mean(Peak_Thrust)
    std_maxThurst(N) = std(Peak_Thrust)

    avg_Impules(N) = trapz(tot_data)
    
    %Calibration = w_avg_std(dataN,uncertainty,Time);
    %avg_curve(N,:) = Calibration(:,2);




    % Note that averaging should accour before data fitting. Technically
    % either can be done, but the output will be much more smooth if the
    % fit is applied to an average
    %% Data Fitting
    data = tot_data;

    [maxi,max_index] = max(data);

    T1 = Time(Time <= Time(max_index));
    T2 = Time(Time > Time(max_index));
    D1 = data((Time <= Time(max_index)));
    D2 = data((Time > Time(max_index)));
    P1 = polyfit(T1,D1,5);
    P2 = polyfit(T2,D2,13);
    fit1 = polyval(P1,T1);
    fit2 = polyval(P2,T2);
    
    tot_fit = zeros(1,501);
    for i = 1:length(fit1)
        tot_fit(i) = fit1(i);
    end
    for i = 1:length(fit2)
        tot_fit(i+length(fit1)) = fit2(i);
    end
    model(N,:) = tot_fit;
    % figure(N);
    % plot(Time,tot_fit)
    % hold on; 
    % plot(Time,data)



    %% Sample onto the standard output array format
    
    


% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
    %% Convert to table for output
    % It is very important that the data is 501 elements long corresponding
    % to 0-0.5 seconds of time at this point!!!
    ThrustCurves(:, N) = model(N,:);
    % Header naming convention of <bottle size (in ml)>_<water volume (in ml)>
    ThrustCurvesNames{N} = [bottleSize, '_', num2str(waterSize)];
end
ThrustCurves = array2table(ThrustCurves);
ThrustCurves.Properties.VariableNames = ThrustCurvesNames;
end
