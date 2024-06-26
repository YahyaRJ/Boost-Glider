function [dSdt] = BoostAscent_odefun(t,S,consts,thrustVec,Time)
% BOOSTASCENT_odefun
% The goal of this function is to output the derivatives of the current
% state (held in the vector 'S'). These derivatives are calculated using
% fundamental physics, and then are packed together in "dSdt'. 
% 
% ODE45 will pass in time information to this function (since the physics
% are time dependentant), and then will numerically integrate the 
% derivative in the output by something like:
%       S(i+1) = S(i) + dSdt*dt 
% where ODE45 is picking dt.
%
% We need information about the position, velocity, and mass, thus these
% are all of the quantities held in the state vector. Also passed in is a
% constants vector that holds information that the will be needed to
% calculate quanties of interest in this funciton. Finally, there is also
% the trust vector and time vector associated with the thrust vector so
% that we can calculate thrust at any given time using interpolation.

%% Unpack the state vector
Vx = S(1); % inertial velocity in x-direction [m/s]
Vy = S(2); % inertial velocity in y-direction [m/s]
Vz = S(3); % inertial velocity in z-direction [m/s]
x  = S(4); % position in x (inertial) [m]
y  = S(5); % position in y (inertial) [m]
z  = S(6); % position in z (inertial) [m]
m  = S(7); % current total mass [kg]

%% Unpack Constants Vector
% Basic Properties
g       = consts(1); % Accelatation due to gravity [m/s^2]
rho_w   = consts(2); % Density of water [kg/m^3]
rho_a   = consts(3); % Density of air [kg/m^3]
mu_k    = consts(4); % Launch rail coefficient of dynamic friction []
% Vehicle info
A_exit  = consts(5); % Area of bottle outlet [m^2]
C_D     = consts(6); % C_D of the vehicle (assume zero lift) []
S_ref   = consts(7); % Wing reference area [m^2]
m_empty = consts(8); % Weight of the rocket with no water [kg]
% Wind
Wx      = consts(9); % Inertial wind velocity in x [m/s] 
Wy      = consts(10); % Inertial wind velocity in x [m/s] 
% Launch Direction
eliv    = consts(11); % Launch elevation [degrees]
azim    = consts(12); % Launch Azimuth [degrees], measured CW from north when looking down on the map; also known as compass heading

% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
%% Wind Triangle
Vax = Vx - Wx;% x air-relative velocity (inertial frame)

Vay = Vy - Wy;% y air-relative velocity (inertial frame)

Va = sqrt(Vax^2 + Vay^2); % Airspeed (magnitude, body frame)

%% Set velocity direction
% Constant if still on the rails
f_rails = 0; % friction due to the rails, preallocate to zero
if sqrt(x^2+y^2+z^2) <= 0.5 % still on rails
    % We need to add 90 degrees to elevation as it is normally measured
    % from the horizon but the spherical to cartiesian coordinate
    % conversion measures from the positive-z (strait down)
    hBod_x = sind(eliv+90)*cosd(azim); % heading of the body (inertial pointing direction of the body)
    hBod_y = sind(eliv+90)*sind(azim);
    hBod_z = cosd(eliv+90);
    f_rails = mu_k*m*g; % making friction non-zero only if on the rails, assuming the force on the rails is the full weight (this is a strange assumption but attempts to account for extra forces during launch)
else %free flight, velocity is into the relative wind
    hVrel = 1/Va*[Vax; Vay; Vz]; % relative velocity (body frame) unit vector
    hBod_x = hVrel(1);
    hBod_y = hVrel(2);
    hBod_z = hVrel(3);
end

%% Interpolate from the thrust curve
% We want to do this so that we allow ode 45 to choose its own time step size
if t < 0.5
    T = interp1(Time, thrustVec, t);
else
    T = 0;
end

%% Calculate total drag
dynamic_pressure = 0.5 .* rho_a .* (Va.^2); %checks out
D = dynamic_pressure .* C_D .* S_ref; % checks out

%% Sum the forces
% Assume that all forces exept gravity act in (or against) the direction
% that the body is pointing. This means that the force in any component
% direction is the sum of the forces multiplied by the component of the
% unit vector in line with the body pointing (hBod). We then assume that
% gravity acts only in -z.

% Sum of the forces
Fx = hBod_x .* (T-D);
Fy = hBod_y .* (T-D);
Fz = (hBod_z .* (T-D)) +  (m.*g);

%% Calculate the derivatives we need to pass out
% Change in postition = Acceleration 
dVdt_x = Fx ./ m;
dVdt_y = Fy ./ m;
dVdt_z = Fz ./ m;

% Now we do some error checking:
% If we are not yet producing enough thrust to get a positive component of
% acceleration while very near the base of the the launcher, this means
% that the rocket is trying to slide backwards down the rails. There is
% backing plate that will stop this from happening, so if this backwards
% motion is detected by the following test, set the acceleration to zero
if sqrt(x^2+y^2+z^2) <= 0.05 && dVdt_z > 0% still on rails
    dVdt_x = 0;
    dVdt_y = 0;
    dVdt_z = 0;
end

% Change in postition = Velocity
dxdt = Vx;
dydt = Vy;
dzdt = Vz;

% Mass flow rate
if m > m_empty % if water is not  yet exausted as measured by weight
    mDot = -sqrt(rho_w .* A_exit .* abs(T));  %stupid, had wrong sign
    %mDot = 1;      
else % ignore any mass change from expulsed air
    mDot = 0;
end
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////

%% Pass out the derivatives in time for ODE45 to intake
% This needs to be in the same 
dSdt = [dVdt_x; dVdt_y; dVdt_z; dxdt; dydt; dzdt; mDot];
end


%{


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

%}