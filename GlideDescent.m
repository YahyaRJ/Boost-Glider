function [GlideData] = GlideDescent(LD, apogee, Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count)

% [GlideRange] = GlideDescent(LD_mod1, apogee, Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count);

%% GlideDescent Summary:
% This funciton is meant to take in and process L/D data and boost apogee
% to find a best glide range, along with wind data to find actual glide 
% range. These assume no wind and an instantanious change to best glide
% velocity and trim. Note that while there are three oswalds models and
% thus three L/D curves, this function only takes one in. This is because
% at this point, a final oswalds model should be chosen, and the L/D curve
% for that model should be the only one passed in to this funciton.
%
% Also note that the heading for best glide should be the apopgee heading
% (as there is no wind, so no reason for the heading to change), while for
% actual glide we will assume the aircraft immediatly turns directly into
% the true wind (inertial wind direction) and holds that heading

%% Outputs:
% GlideData:
%   A table of scalar values of the best glide range from apogee in [m] and
%   key flight parameters for that best glide range including glide angle
%   (theta), CL for best glide (CL_LDmax), velocity at best glide 
%   (V_LDmax), sink velocity (Vsink), and the AoA required for best glide
%   (AoA_LDmax) for each case input.

%% Preallocate variables of interest
bestGlide = zeros(Count,1); % Best glide range from apogee [m]
LDmax = zeros(Count,1); % L/D max value used []
theta = zeros(Count,1); % Glide angle [degrees]
CL_LDmax = zeros(Count,1); % CL required for best glide []
V_LDmax = zeros(Count,1); % Velocity required for best glide [m/s]
Vsink = zeros(Count,1); % Sink rate at best glide [m/s]
AoA_LDmax = zeros(Count,1); % AoA required for best glide [degrees]
WingLoading = zeros(Count,1); % Wing loading at best glide [N/m^2]

%% Loop through different configurations
for n = 1:Count

    % /////////////////////////////////////////////////////////////////////////
    % MODIFY THIS SECTION
    % /////////////////////////////////////////////////////////////////////////
    %% Find L/D max
    [M,I] = max(LD{n,:});  %get the index for C_L_L/D_max
    LDmax(n) = M;

    %% Find best glide range
    bestGlide(n) = LDmax(n) .* apogee(n);

    %% GLide angle
    theta(n) = atand(1 ./ LDmax(n)); %in degrees

    %% Find Best Glide CL and Velocity
    CL_LDmax(n) = WingLiftCurve{n,I};
    %subtract water and compressed air weight(air not even there?!)
    weight = Weight_Data.Wo(n) - Weight_Data.W_water(n);
    V_LDmax(n) = sqrt((2.*weight) ./ (ATMOS.rho(n) .* CL_LDmax(n) .* Design_Input.Sref_w(n)));   
    Vsink(n) = V_LDmax(n) .* sind(theta(n));
    angles_of_attack = -5:12;
    alpha = angles_of_attack(I); %angle of attack at which LD is max
    AoA_LDmax(n) = (CL_LDmax(n) ./ alpha) + WingLiftModel.AoA_0(n);

    %% Find Wing Loading
    WingLoading(n) = weight ./ Design_Input.Sref_w(n);
    % /////////////////////////////////////////////////////////////////////////
    % END OF SECTION TO MODIFY
    % /////////////////////////////////////////////////////////////////////////
end

%% Convert to tables for output
GlideData = table(bestGlide, LDmax, theta, CL_LDmax, V_LDmax, Vsink, AoA_LDmax, WingLoading);

end

