%{
Author(s): Ghulam Yahya Rajabi and Zachary Purdey
Assignment title: Project 2 Group Code
Purpose: Determine the effect each parameter(initial water volume, drag
coefficient, initial air pressure, and the launch angle) has on the
furthest distance that the rocket can travel and by varying these
parameters, determine the combination of parameters which would allow to
hit a target 85meters away
Creation date: 12/04/2023
Revisions: Please see below

%}


clear;
clc;
close all;


% get the constants
const = getConst();

%load the verification
load('project2verification.mat');



%revision #1: Date: 12/05 Yahya calculated the effect of each parameter
%% Drag Coefficient Parameters

%drag can be from 0.3 to 0.5
drag_parameters = 0.3:0.01:0.5;

%store original drag
drag_original = const.C_D;

%variable to store the max distances
drag_distance = zeros(size(drag_parameters));

for i = 1:length(drag_parameters)

    % vary the drag coefficient
    const.C_D = drag_parameters(i);

    % intitial conditions
    initial_conditions = [const.x_0, const.v_x_0, const.z_0, ...
    const.v_z_0, const.m_i_r, const.V_i_air, const.m_i_air]; 

    % ode45
    [t, X] = ode45(@(t, X) RocketODE(t, X, const), const.integrationTime, ...
    initial_conditions);

    %store max distance
    drag_distance(i) = max(X(X(:,3)>0 ,1));
end

% Find the index at which distance is maximum
[maxDistance_drag, max_index_drag] = max(drag_distance);
% Find the optimal value of this parameter
drag_max = drag_parameters(max_index_drag);

%plot it
figure(1);
subplot(4, 2, 1);
plot(drag_parameters, drag_distance);
title('Drag Coefficient Variation');
xlabel('Drag Coefficient');
ylabel('Distance (m)');
xline(0.3, 'r--', 'LineWidth', 1);
xline(0.4, 'b--', 'LineWidth', 1);

%reset the drag coefficient
const.C_D =  drag_original;

%% Angle of Attack Parameters

%angle can vary from 0 to 90 degrees, 
angle_parameters = 0:1:90;

%store original angle
angle_original = const.theta_i;

%variable to store the max distances
angle_distance = zeros(size(angle_parameters));

for i = 1:length(angle_parameters)

    % vary the angle
    const.theta_i = angle_parameters(i);

    % intitial conditions
    initial_conditions = [const.x_0, const.v_x_0, const.z_0, ...
    const.v_z_0, const.m_i_r, const.V_i_air, const.m_i_air]; 

    % ode45
    [t, X] = ode45(@(t, X) RocketODE(t, X, const), const.integrationTime, ...
    initial_conditions);

    %store max distance
    angle_distance(i) = max(X(X(:,3)>0 ,1));
end

% Find the index at which distance is maximum
[maxDistance_angle, max_index_angle] = max(angle_distance);
% Find the optimal value of this parameter
angle_max = angle_parameters(max_index_angle);

%plot it
subplot(4, 2, 2);
plot(angle_parameters, angle_distance);
title('Launch Angle Variation');
xlabel('Angle (degrees)');
ylabel('Distance (m)');
xline(angle_max, 'r--', 'LineWidth', 1);

%reset the angle
const.theta_i =  angle_original;


%% Pressure Parameters

%gauge pressure from 0 to 150 psi, converted to Pa
pressure_parameters = 0:1:150; %psi
pressure_parameters = 6894.76 .* pressure_parameters; %Pa

%store original values
gauge_pressure_original = const.p_00;
pressure_original = const.p_0;
massrocket_original = const.m_i_r;
massair_original = const.m_i_air;
endpressure_original = const.P_air_end;

%variable to store the max distances
pressure_distance = zeros(size(pressure_parameters));

for i = 1:length(pressure_parameters)

    % vary the pressure
    const.p_00 = pressure_parameters(i); % gauge pressure
    const.p_0 = const.p_00 + const.p_a; %Pa: intitial pressure of air in bottle
    %initial mass of bottle
    const.m_i_r = const.m_B + (const.rho_w .* const.V_i_water) + ...
    ((const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0));
    %initial mass of air
    const.m_i_air = (const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0);
    %pressure of air at end
    const.P_air_end = const.p_0*((const.V_i_air/const.V_B)^const.gamma);


    % intitial conditions
    initial_conditions = [const.x_0, const.v_x_0, const.z_0, ...
    const.v_z_0, const.m_i_r, const.V_i_air, const.m_i_air]; 

    % ode45
    [t, X] = ode45(@(t, X) RocketODE(t, X, const), const.integrationTime, ...
    initial_conditions);

    %store max distance
    pressure_distance(i) = max(X(X(:,3)>0 ,1));
end

% Find the index at which distance is maximum
[maxDistance_pressure, max_index_pressure] = max(pressure_distance);
% Find the optimal value of this parameter
pressure_max = pressure_parameters(max_index_pressure);

%plot it
subplot(4, 2, 3);
plot(pressure_parameters, pressure_distance);
title('Gauge Pressure Variation');
xlabel('Pressure (Pa)');
ylabel('Distance (m)');
xline(150*6894.76, 'r--', 'LineWidth', 1);
xline(100*6894.76, 'g--', 'LineWidth', 1);

%reset the original values
const.p_00 = gauge_pressure_original;
const.p_0 = pressure_original;
const.m_i_r = massrocket_original;
const.m_i_air = massair_original;
const.P_air_end = endpressure_original;


%% Water Volume Parameters

%volume of water from 0 to 0.002m^3
volume_parameters = 0:0.0001:0.002; 

%store original volume
volumewater_original = const.V_i_water;
volumeair_original = const.V_i_air;
endTemp_original = const.T_air_end;

%variable to store the max distances
volume_distance = zeros(size(volume_parameters));

for i = 1:length(volume_parameters)

    % vary the volume
    const.V_i_water = volume_parameters(i);
    const.V_i_air = const.V_B - const.V_i_water; % initial volume of air
    %initial mass of bottle
    const.m_i_r = const.m_B + (const.rho_w .* const.V_i_water) + ...
    ((const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0));
    %initial mass of air
    const.m_i_air = (const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0);
    %Pressure of air at the end when no water left
    const.P_air_end = const.p_0*((const.V_i_air/const.V_B)^const.gamma);
    const.T_air_end = const.T_0*((const.V_i_air/const.V_B)^(const.gamma - 1));

    % intitial conditions
    initial_conditions = [const.x_0, const.v_x_0, const.z_0, ...
    const.v_z_0, const.m_i_r, const.V_i_air, const.m_i_air]; 

    % ode45
    [t, X] = ode45(@(t, X) RocketODE(t, X, const), const.integrationTime, ...
    initial_conditions);

    %store max distance
    volume_distance(i) = max(X(X(:,3)>0 ,1));
end

% Find the index at which distance is maximum
[maxDistance_volume, max_index_volume] = max(volume_distance);
% Find the optimal value of this parameter
volume_max = volume_parameters(max_index_volume);

%plot it
subplot(4, 2, 4);
plot(volume_parameters, volume_distance);
title('Water Volume Variation');
xlabel('Volume (m^3)');
ylabel('Distance (m)');
xline(volume_max, 'r--', 'LineWidth', 1);

%reset the original values
const.V_i_water = volumewater_original;
const.V_i_air = volumeair_original;
const.T_air_end = endTemp_original;
const.m_i_r = massrocket_original;
const.m_i_air = massair_original;
const.P_air_end = endpressure_original;


% revision 2: Date 12/06/2023   Zach calculated the parameters to hit
% target distance of 85m

%% Method 1: 

%assume a factor of safety of 1.5 and use 100 psi for the gauge pressure
%take the median of the 0.3-0.5 coefficient of drags
%since, initial volume of water is an easily manipulable parameter, use the
%optimum value
%get the angle at which target is hit

%change the constants
const.C_D = 0.4;
const.p_00 = 100 * 6894.76; %100 psi converted to Pa
const.V_i_water = volume_max;
%to change them, below must also change
const.p_0 = const.p_00 + const.p_a;
const.V_i_air = const.V_B - const.V_i_water;
const.m_i_r = const.m_B + (const.rho_w .* const.V_i_water) + ...
((const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0));
const.m_i_air = (const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0);
const.P_air_end = const.p_0*((const.V_i_air/const.V_B)^const.gamma);
const.T_air_end = const.T_0*((const.V_i_air/const.V_B)^(const.gamma - 1));

% intitial conditions
initial_conditions = [const.x_0, const.v_x_0, const.z_0, ...
const.v_z_0, const.m_i_r, const.V_i_air, const.m_i_air];

%angle can vary from 0 to 90 degrees, 
new_angle_parameters = 0:0.1:90;

%variable to store the distances
new_angle_distance = zeros(size(angle_parameters));

for i = 1:length(new_angle_parameters)

    % vary the angle
    const.theta_i = new_angle_parameters(i); 

    % ode45
    [t, X] = ode45(@(t, X) RocketODE(t, X, const), const.integrationTime, ...
    initial_conditions);

    %store max distance
    new_angle_distance(i) = max(X(X(:,3)>-0.7 ,1));
end

%our target angle
launchAngleForTargetDistance = interp1(new_angle_distance(1:200), new_angle_parameters(1:200), 85, 'linear', 'extrap');

%set the theta as that angle
const.theta_i = launchAngleForTargetDistance;

%run the ODE again
[t, X] = ode45(@(t, X) RocketODE(t, X, const), const.integrationTime, ...
    initial_conditions);

% Plot Height vs. Distance
d = X(:,1); % distance (x position)
h = X(:,3); % height (z position)
figure(2)
plot(d(h>-1),h(h>-1),'LineWidth',2);
title('Rocket Trajectory');
ylabel('Height (m)');
xlabel('Distance (m)');
yline(0);
xline(84, 'r--', 'LineWidth', 1);
text(87, 5, '85 \pm 1', 'HorizontalAlignment', 'left', 'Color', 'red');
xline(86, 'r--', 'LineWidth', 1);
xlim([0,100]);
ylim([0,inf]);
text(50, 1, '\theta = 12.2^\circ', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(50, 4, 'P_{gauge} = 100.0 psi', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(50, 3, ['V = ' sprintf('%.1f', volume_max*1000) ' L'], 'HorizontalAlignment', 'center', 'FontSize', 12);
text(50, 2, 'C_D = 0.4', 'HorizontalAlignment', 'center', 'FontSize', 12);

% revision 4: Date 12/10/2023   Yahya calculated the maximum distance
% possible

%% Get the Maximum Distance
%find the optimal values for each parameter and run the ODE with it to see
%how far the rocket can reach

%change the constants
const.theta_i = angle_max; 
const.C_D = drag_max;
const.p_00 = pressure_max; 
const.V_i_water = volume_max;
%to change them, below must also change
const.p_0 = const.p_00 + const.p_a;
const.V_i_air = const.V_B - const.V_i_water;
const.m_i_r = const.m_B + (const.rho_w .* const.V_i_water) + ...
((const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0));
const.m_i_air = (const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0);
const.P_air_end = const.p_0*((const.V_i_air/const.V_B)^const.gamma);
const.T_air_end = const.T_0*((const.V_i_air/const.V_B)^(const.gamma - 1));

% intitial conditions
initial_conditions = [const.x_0, const.v_x_0, const.z_0, ...
const.v_z_0, const.m_i_r, const.V_i_air, const.m_i_air];

%run the ODE again
[t, X] = ode45(@(t, X) RocketODE(t, X, const), [0,200], ...
    initial_conditions);

% Plot Height vs. Distance
d = X(:,1); % distance (x position)
h = X(:,3); % height (z position)
figure(3)
plot(d(h>-20),h(h>-20),'LineWidth',2);
title('Rocket Trajectory');
ylabel('Height (m)');
xlabel('Distance (m)');
yline(0);
xline(177.45, 'r--', 'LineWidth', 2);
ylim([0, inf]);
xlim([0, 200]);
text(175, 62, '177.45 \pm 1', 'HorizontalAlignment', 'right', 'Color', 'red');
text(100, 50, ['P_{gauge} = ' sprintf('%.1f', pressure_max/6894.76) ' psi'], 'HorizontalAlignment', 'center', 'FontSize', 12);
text(100, 40, ['V = ' sprintf('%.1f', volume_max*1000) ' L'], 'HorizontalAlignment', 'center', 'FontSize', 12);
text(100, 30, ['C_D = ' sprintf('%.1f', drag_max)], 'HorizontalAlignment', 'center', 'FontSize', 12);
text(100, 20, ['\theta = ' sprintf('%.1f', angle_max) '^\circ'], 'HorizontalAlignment', 'center', 'FontSize', 12);


%% part 1

%{

% intitial conditions
initial_conditions = [const.x_0, const.v_x_0, const.z_0, ...
    const.v_z_0, const.m_i_r, const.V_i_air, const.m_i_air]; 

% ode45
[t, X] = ode45(@(t, X) RocketODE(t, X, const), const.integrationTime, ...
    initial_conditions);

thrust = thrustFunc(t, X(:,6), X(:,7), const);

%% plots
d = X(:,1); % distance (x position)
h = X(:,3); % height (z position)

% Plot Height vs. Distance
figure(1)
subplot(2, 1, 1);
plot(d(h>-2),h(h>-2),'LineWidth',2);
hold on;
plot(verification.distance, verification.height, '--', 'LineWidth', 2);
title('Rocket Trajectory');
ylabel('Height (m)');
xlabel('Distance (m)');
xline(2.9341,'--','Color','b');
xline(4.3053,'--', 'Color', 'r');
ylim([0,inf]);
subplot(2, 1, 2);
plot(t,thrust, 'LineWidth', 2);
hold on;
plot(verification.time, verification.thrust, '--', 'LineWidth', 2);
title('Thrust vs Time');
ylabel('Thrust (N)');
xlabel('T (s)');
xline(0.1837,'--','Color','b');
xline(0.2159,'--', 'Color', 'r');
print(gcf,'verification.png','-dpng','-r800');

%}




%% Functions

% Constants
function const = getConst()
    const.g = 9.807;  % m/s^2: acceleration due to gravity
    const.c_dis = 0.8; % discharge coefficient
    const.rho_air = 0.961; % kg/m3: ambient air density
    const.V_B = 0.002; % m3: volume of empty bottle
    const.p_a = 83426.563088; % Pa: atmospheric pressure
    const.gamma = 1.4; % ratio of specific heats for air
    const.rho_w = 1000; %kg/m^3: water density
    const.d_e = 2.1; %cm: diameter of the throat
    const.a_e = pi/4*(const.d_e/100)^2; % m^2 cross sectional area of throat
    const.d_B = 10.5; %cm: diameter of the bottle
    const.a_B = pi/4*(const.d_B/100)^2; % m^2 cross sectional area of bottle
    const.R_air = 287; %J/(kg.K): specific gas constant of air
    const.m_B = 0.15; %kg: mass of empty 2-liter bottle with cone and fins
    const.C_D = 0.48; % drag coefficient
    const.p_00 = 358527.37856; %Pa, given gauge pressure
    const.p_0 = const.p_00 + const.p_a; %Pa: intitial pressure of air in bottle
    const.V_i_water = 0.00095; %m^3: initial volume of water inside bottle
    const.T_0 = 300; %K: initial temperature of air
    const.v_0 = 0.0; %m/s: initial velocity of rocket
    const.v_x_0 = 0; %m/s initial in x direction
    const.v_z_0 = 0; %m/s initial in z direction
    const.theta_i = 42; %initial angle of rocket, degrees
    const.x_0 = 0.0; % m: initial horizontal distance
    const.z_0 = 0.25; % m: initial vertical distance
    const.l_s = 0.5; %m: length of test stand
    const.integrationTime = [0 5]; %integration time
    const.V_i_air = const.V_B - const.V_i_water; % initial volume of air
    %initial mass of bottle
    const.m_i_r = const.m_B + (const.rho_w .* const.V_i_water) + ...
    ((const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0));
    %initial mass of bottle
    const.m_i_air = (const.p_0 .* const.V_i_air)./(const.R_air .* const.T_0);
    %Pressure of air at the end when no water left
    const.P_air_end = const.p_0*((const.V_i_air/const.V_B)^const.gamma);
    const.T_air_end = const.T_0*((const.V_i_air/const.V_B)^(const.gamma - 1));
end

function ddt = RocketODE(t, X, const)

    %calculate air pressure inside to differentiate between Stage 2 or 3
    P_air = ((X(7) ./ const.m_i_air).^const.gamma) .* const.P_air_end;
    
    if X(6) < const.V_B 
        %Stage 1
        
        %calculate pressure
        P = ((const.V_i_air/X(6))^const.gamma)*const.p_0;

        %no change in mass of air
        ddt(7) = 0;

        % change in mass of rocket
        v_e = sqrt( (2*(P - const.p_a)) ./ (const.rho_w));
        ddt(5) = -const.c_dis * const.rho_w * const.a_e * v_e;

        %Equation 10:
        ddt(6) = const.c_dis .* const.a_e .* sqrt((2./const.rho_w).* ...
            (const.p_0.*((const.V_i_air/X(6))^const.gamma)-const.p_a));

        %thrust
        F = 2*const.c_dis*(P-const.p_a)*const.a_e;

    elseif P_air > const.p_a
        %Stage 2
        
        %no change in volume of air
        ddt(6) = 0;

        %calculate density and temperature of air inside
        rho_air_inside = X(7) / const.V_B;
        T_air_inside = P_air ./ (const.R_air .* rho_air_inside);
        
        %ciritical pressure
        P_critical = P_air*((2/(const.gamma+1))^(const.gamma/(const.gamma-1)));

        %flow choked?
        if P_critical > const.p_a
            %choked

            %exit pressure, temperature, density and velocity
            P_exit = P_critical;
            T_exit = (2/(const.gamma+1))*T_air_inside;
            rho_exit = P_exit/(const.R_air*T_exit);
            v_exit = sqrt(const.gamma .* const.R_air .* T_exit);

        else
            %not choked

            %exit pressure, temperature, density and velocity
            P_exit = const.p_a;
            %exit Mach number
            %Mach_exit = sqrt((2/(const.gamma - 1)) * ...
            % (exp((const.gamma - 1)/const.gamma * log(P_air/const.p_a)) - 1));
            Mach_exit = sqrt((2/(const.gamma - 1)) * ((P_air/const.p_a)^(0.4/1.4) - 1));
            T_exit = T_air_inside / (1+((const.gamma - 1)/2)*(Mach_exit^2));
            rho_exit = const.p_a/(const.R_air*T_exit);
            v_exit = Mach_exit .* sqrt(const.gamma .* const.R_air .* T_exit);

        end
        
        %air mass flow
        m_dot_air = const.c_dis * rho_exit * const.a_e * v_exit;

        %change in mass of air in the rocket, equation 21
        ddt(7) = -m_dot_air;

        %change in mass of rocket
        ddt(5) = ddt(7);

        %thrust 
        F = (m_dot_air * v_exit) + (P_exit - const.p_a)*const.a_e;

    else
        %Stage 3

        %no change in mass of bottle
        ddt(5) = 0;
        %no change in volume of air
        ddt(6) = 0;
        %no change in mass of air
        ddt(7) = 0;
        %no thrust generated
        F = 0;
    end
    
    %figure out the heading
    if sqrt(X(1)^2 + (X(3) - const.z_0)^2) < const.l_s
        % On test stand
        heading_x = cosd(const.theta_i); 
        heading_z = sind(const.theta_i);
    else
        % Not on test stand
        heading_x = X(2) / sqrt(X(2)^2 + X(4)^2);
        heading_z = X(4) / sqrt(X(2)^2 + X(4)^2);
    end
    %h_hat = heading_x + heading_z;

    %drag
    v_h = sqrt(X(2)^2 + X(4)^2);
    D = (const.rho_air/2)*(v_h^2) * const.C_D * const.a_B;

    % velocity in the x direction
    ddt(1) = X(2);
    % velocity in the z direction
    ddt(3) = X(4);
    % acceleration in the x direction
    ddt(2) = ((F-D)*heading_x)/X(5);
    % acceleration in the z direction
    ddt(4) = (((F-D)*heading_z)/X(5))-const.g;
    
    % transpose for ODE45
    ddt = ddt';
    
end

function thrust = thrustFunc(t, V_air, m_air, const)
    for i = 1:length(t)
        %pressure of air inside rocket assuming stage 2 or 3
        P_air_inside = ((m_air(i) / const.m_i_air)^const.gamma)*const.P_air_end;
        if V_air(i) < const.V_B
        %phase 1
            %pressure inside
            P = ((const.V_i_air/V_air(i))^const.gamma)*const.p_0;
            % thrust
            thrust(i) = 2* const.c_dis*(P-const.p_a)*const.a_e;
        elseif P_air_inside > const.p_a
            %air density
            rho = m_air(i)/const.V_B;
            % temperature
            T = P_air_inside/(rho*const.R_air);
            % critical pressure
            P_critical = P_air_inside*((2/(const.gamma+1))^ ...
                (const.gamma/(const.gamma-1)));
            %if choked
            if P_critical > const.p_a
                % exit P, T, rho, v
                P_exit = P_critical;
                T_exit = (2/(const.gamma+1))*T;
                rho_exit = P_exit/(const.R_air*T_exit);
                v_exit = sqrt(const.gamma*const.R_air*T_exit);
            %not choked
            else
                % exit P, T, rho, v
                P_exit = const.p_a;
                Mach_exit = sqrt((2/(const.gamma - 1)) * ((P_air_inside/const.p_a)^(0.4/1.4) - 1));
                T_exit = T / (1+((const.gamma - 1)/2)*(Mach_exit^2));
                rho_exit = const.p_a/(const.R_air*T_exit);
                v_exit = Mach_exit .* sqrt(const.gamma .* const.R_air .* T_exit);
            end
            %mass flow
            m_dot_air = const.c_dis * rho_exit * const.a_e * v_exit;

            %thrust 
            thrust(i) = (m_dot_air * v_exit) + (P_exit - const.p_a)*const.a_e;

        else
        %phase 3
            thrust(i) = 0;
        end
            
    end
end

