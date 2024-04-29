clear
clc
close all

%load data
data = load('boost2');
T = data(:,1);
% We have to convert pressure P from mBar to Pa,and input to atmospalt
P = data(:,2);
h = atmospalt(100*P); %height in m
init_alt = h(1);
% Subtract the initial altitude out so that relative altitude is displayed
H = h - init_alt; %Height diff in m

%630-936

figure();
plot(T(630:936),H(630:936));
xlabel('Time (s)');
ylabel('Apogee (m)');
title('Boost Ascent Test #2');


%load data
data2 = load('boost3');
T2 = data2(:,1);
% We have to convert pressure P from mBar to Pa,and input to atmospalt
P2 = data2(:,2);
h2 = atmospalt(100*P2); %height in m
init_alt2 = h2(1);
% Subtract the initial altitude out so that relative altitude is displayed
H2 = h2 - init_alt2; %Height diff in m

%523-915

figure();
plot(T2(523:915),H2(523:915));
xlabel('Time (s)');
ylabel('Apogee (m)');
title('Boost Ascent Test #3');