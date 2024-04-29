clear
clc
close all

%load data
data = load('boost3');
T = data(:,2)
% We have to convert pressure P from mBar to Pa,and input to atmospalt
P = data(:,2);
h = atmospalt(100*P); %height in m
init_alt = h(1);
% If this is the first loop save the initial altitude so that everything else is relative
if i == 1
init_alt = h;
end
% Subtract the initial altitude out so that relative altitude is displayed
H = h - init_alt; %Height diff in m

figure();
plot(H);