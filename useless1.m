clear
clearvars
clc
close all

wingspan_values = 0.01:0.01:1.0;
wing_planform_area_values = 0.01:0.01:1.0;

% Generate combinations
combinations = combvec(wingspan_values, wing_planform_area_values);

combinations = combinations';

% Calculate aspect ratio for each combination
aspect_ratio = combinations(:,1).^2 ./ combinations(:,2);

% Add aspect ratio as the third column
combinations = [combinations, aspect_ratio];

%disp(combinations);

a = [1,2,3];
b= [4,5,6];
comb = combvec(a,b);
comb = comb';
