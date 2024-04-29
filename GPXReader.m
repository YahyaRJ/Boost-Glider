%% GPX File Reader
% ASEN 2804
% Author: Bryce Pfuetze
% Date: Initiated - 20 Dec. 2021
% Last modified - 3 April 2024

% Step 1: input your file name
% Step 1: run this script for the first time, it will plot your whole run
% Step 2: mark the start and end of your trial (the Elevation vs Index 
%         plot helps). Set 'beginning' and 'ending' equal to these indices
%         and uncomment them.
% Step 3: Set 'cropped' equal to 1
% Step 4: Run the script again on your cropped data

clc; clear; close all

% Input file name here!
%fileName = 'glide1.gpx';
fileName = 'glide2.gpx';

%beginning = 120;
%ending = 131;

beginning = 107;
ending = 116;

cropped = 1; % uncropped

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load GPX data into a structure
trial = readstruct(fileName, 'FileType','xml','StructNodeName','trkseg');

% Pull values from the structure
latitude = [trial.trkpt.latAttribute]; % deg
longitude = [trial.trkpt.lonAttribute]; % deg
elevation = [trial.trkpt.ele]; % meters

t = length(trial.trkpt); % 1 Hz sampling
time = 1:t; % seconds

% Convert Spherical to Cartesian
theta = zeros(1,t);
phi = zeros(1,t);

earthRad = 6369339; % meters, in Boulder (40 deg North)
realElevation = elevation + earthRad;
xcoord = realElevation.*cosd(latitude).*cosd(longitude);
ycoord = realElevation.*cosd(latitude).*sind(longitude);
zcoord = realElevation.*sind(latitude);

if cropped == 1
    xcoord = xcoord(beginning:ending);
    ycoord = ycoord(beginning:ending);
    zcoord = zcoord(beginning:ending);
    elevation = elevation(beginning:ending);
    latitude = latitude(beginning:ending);
    longitude = longitude(beginning:ending);
    t = ending - beginning+1;
    time = 1:t;
else
    fprintf("The data has not been cropped")
end


% You can also use the Automated Driving Toolbox.
% This takes into account mapping projections to make everything line up
% as we would expect. For our purposes, assuming these trials are over a
% short path, basic spherical conversions are sufficient. You can
% uncomment the following two lines and line 66 to see a comparison of
% the two.
% origin = [latitude(1), longitude(1), elevation(1)];
% [xcoord1,ycoord1, zcoord1] = latlon2local(latitude,longitude,elevation,origin);


% Calculate displacements
horizontalDisp = zeros(1,t);
verticalDisp = zeros(1,t);
totalDisp = zeros(1,t);

for i=1:(t-1)
    verticalDisp(i) = elevation(i + 1) - elevation(i);
    horizontalDisp(i) = sqrt((xcoord(i+1) - xcoord(i))^2 + (ycoord(i+1) - ycoord(i))^2);
    totalDisp(i) = sqrt((xcoord(i+1) - xcoord(i))^2 + (ycoord(i+1) - ycoord(i))^2 + (zcoord(i+1) - zcoord(i))^2);
end

% Use displacement to find ground track length, speed, and gradient (%)
distance = zeros(1,t);
for i=1:(t-1)
    distance(i+1) = horizontalDisp(i+1) + distance(i);
end

speed = totalDisp / 1; % m/s with 1 Hz sampling
grade = verticalDisp./horizontalDisp .*100; % gradient as a percent
elevMin = min(elevation);

figure
area(cumsum(horizontalDisp)-horizontalDisp(1),elevation,elevMin);
ylabel('Elevation (meters)')
xlabel('Horizontal Displacement [m]')
title('Elevation');

figure
area(1:length(elevation),elevation,elevMin);
ylabel('Elevation (meters)')
xlabel('Index')
title('Elevation vs Index');

figure
plot3(xcoord - xcoord(1),ycoord-ycoord(1), elevation,'b');
%plot3(xcoord - xcoord(1),ycoord-ycoord(1), elevation,'b', xcoord1 - xcoord1(1),ycoord1-ycoord1(1), elevation);
hold on
plot3(0,0, elevation(1), "og");
plot3(xcoord(end)-xcoord(1),ycoord(end)-ycoord(1), elevation(end), "or");
ax.DataAspectRatio = [1 1 0.5];
xlabel('X position (East-West) [m]')
ylabel('Y Position (North-South) [m]')
zlabel('Elevation [m]')
title('3D Track')
view([8.9 60.7])
hold off

%webmap()
%wmline(latitude,longitude);
if cropped == 1
    horizontalDist = sum(horizontalDisp)
    verticalDist = elevation(1) - elevation(end)
    glideRatio =  horizontalDist/ verticalDist
    glideAngle = rad2deg(atan(verticalDist / horizontalDist))
end