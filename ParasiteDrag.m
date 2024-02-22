function [Parasite_Data,FF_Table] = ...
    ParasiteDrag(Design_Input,Airfoil,WingGeo_Data,ATMOS,Count)
%%  Parasite Drag Summary
% This function performs the Raymer Component Drag Buildup Method to
% determine total and component parasite drag coefficients using the
% aircraft configuration geometry information from the Design Input,
% Airfoil, and WingGeo_Data tables.  Additionally, it leverages the 
% standard atmosphere properties from the ATMOS table. The output table 
% for this function includes the total parasite drag coefficient (CDo), and
% a breakdown of the contribution to CDo for the fuselage (f), wing (w), 
% horizontal stabilizers (h1, h2), vertical stablizers (v1,v2), misc base 
% drag (misc), and leakage and proturbance drag (lp).  This function also
% outputs the total wetted area for the entire aircraft (Swet_total). 
% Finally, the FF_Table consolidates the form factors (FF) which account
% for zero lift pressure drag for each component for evaluation.
% Note that no FF should be < 1.0.  All these calculation are done for each
% configruation in the Design Input spreadsheet.

%% Outputs:
%
% Parasite Data:
%   Table containing total parasite drag and a component breakdown of
%   contributions to that total (columns), for each input case (rows)
%
% FF_Table:
%   Table containing total form factor and a component breakdown of
%   contributions to that total (columns), for each input case (rows)

%% Preallocate variables of interest
CDo = zeros(Count, 1); % Total parasite drag coefficient
CDo_f = zeros(Count, 1); % Fusalage parasite drag coefficient contribution
CDo_w = zeros(Count, 1); % Wing parasite drag coefficient contribution
CDo_h1 = zeros(Count, 1); % Horizontal stabilizer 1 parasite drag coefficient contribution
CDo_h2 = zeros(Count, 1); % Horizontal stabilizer 2 parasite drag coefficient contribution
CDo_v1 = zeros(Count, 1); % Vertical stabilizer 1 parasite drag coefficient contribution
CDo_v2 = zeros(Count, 1); % Vertical stabilizer 2 parasite drag coefficient contribution
CDo_misc = zeros(Count, 1); % Misc. parasite drag coefficient contribution
CDo_lp = zeros(Count, 1); % Leakage and purturbance parasite drag coefficient contribution
Swet_tot = zeros(Count, 1); % Total wetted area [m^2]

FF_f = zeros(Count, 1); % Fusalage form factor
FF_w = zeros(Count, 1); % Wing form factor
FF_h1 = zeros(Count, 1); % Horizontal stabilizer 1 form factor
FF_h2 = zeros(Count, 1); % Horizontal stabilizer 2 form factor
FF_v1 = zeros(Count, 1); % Vertical stabilizer 1 form factor
FF_v2 = zeros(Count, 1); % Vertical stabilizer 2 form factor

%% Loop through different configurations
for n = 1:Count
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
    %% Essentials 
        M = Design_Input.V_o(n)/ATMOS.a(n);
    %% Fuselage Contribution To CDo
    k_surface = 0.405*10^-5; % Chose from slides - chose sheet metal (m)
    Re_f_L = Design_Input.Length_f(n)*Design_Input.V_o(n)/ATMOS.nu(n); %Re for fuselage
    Re_f_cutoff = 38.21*(Design_Input.Length_f(n)/k_surface)^1.053; %Re model using surface roughness
    if Re_f_L > Re_f_cutoff
        Re_f = Re_f_cutoff;
    else
        Re_f = Re_f_L;
    end
    Cf_f = .455/((log10(Re_f))^2.58*(1+.144*M^2)^.65); %Leaving off mach correction
    f_f = Design_Input.Length_f(n)/sqrt(4/pi*Design_Input.Amax_f(n)); %
    FF_f(n) = .9 + 5/(f_f^1.5) + f_f/400; %Fuselage Form Factor
    CDo_f(n) = Cf_f*Design_Input.Q_f(n)*FF_f(n)*Design_Input.Swet_f(n)*18; %Contribution of Fuselage to CDo
    

    %% Wing Contribution to CDo
    Re_w = WingGeo_Data.MAC_w(n)*Design_Input.V_o(n)/ATMOS.nu(n); %Wing Re
    Cf_w = .074/(Re_w^.2); %Wing Flat Plate Coef of Friction for Turbulent Flow
    FF_w(n) = (1+ .6/Airfoil.X_thick_w(n)*Airfoil.Thick_w(n)+100*Airfoil.Thick_w(n)^4)*...
        (1.35*M^.18*cosd(Design_Input.Sweep_w(n))^.28);  %Wing Form Factor
    CDo_w(n) = Cf_w*FF_w(n)*Design_Input.Q_w(n)*Design_Input.Swet_w(n); %Contribution of Wing to CDo

    %% Horizontal Tail #1 Contribution to CDo
    if Design_Input.Swet_h1(n)~=0 % If this component exists:
        Re_h1 = Design_Input.MAC_h1(n)*Design_Input.V_o(n)/ATMOS.nu(n); %Horz Tail Re
        Cf_h1 = .074/(Re_h1^.2); %Flat Plate Coef of Friction for Turbulent Flow
        FF_h1(n) = (1+ .6/Airfoil.X_thick_h1(n)*Airfoil.Thick_h1(n)+100*Airfoil.Thick_h1(n)^4)*...
        (1.35*M^.18*cosd(Design_Input.Sweep_h1(n))^.28); %Horz Tail Form Factor
        CDo_h1(n) = Cf_h1*FF_h1(n)*Design_Input.Q_h1(n)*Design_Input.Swet_h1(n); %Contribution of Horz Tail 1 to CDo 
    end

    %% Horizontal Tail #2 Contribution to CDo
    if Design_Input.Swet_h2(n)~=0 % If this component exists:
        Re_h2 = Design_Input.MAC_h2(n)*Design_Input.V_o(n)/ATMOS.nu(n); %Horz Tail Re
        Cf_h2 = .074/(Re_h2^.2); %Flat Plate Coef of Friction for Turbulent Flow
        FF_h2(n) = (1+ .6/Airfoil.X_thick_h2(n)*Airfoil.Thick_h2(n)+100*Airfoil.Thick_h2(n)^4)*...
        (1.35*M^.18*cosd(Design_Input.Sweep_h2(n))^.28); %Horz Tail Form Factor
        CDo_h2(n) = Cf_h2*FF_h2(n)*Design_Input.Q_h2(n)*Design_Input.Swet_h2(n); %Contribution of Horz Tail 2 to CDo 
    end

    %% Vertical Tail #1 Contribution to CDo
    if Design_Input.Swet_v1(n)~=0 % If this component exists:
        Re_v1 = Design_Input.MAC_v1(n)*Design_Input.V_o(n)/ATMOS.nu(n);%Horz Tail Re
        Cf_v1 = .074/(Re_v1^.2); %Flat Plate Coef of Friction for Turbulent Flow
        FF_v1(n) = (1+ .6/Airfoil.X_thick_v1(n)*Airfoil.Thick_v1(n)+100*Airfoil.Thick_v1(n)^4)*...
        (1.35*M^.18*cosd(Design_Input.Sweep_v1(n))^.28); %Horz Tail Form Factor
        CDo_v1(n) = Cf_v1*FF_v1(n)*Design_Input.Q_v1(n)*Design_Input.Swet_v1(n); %Contribution of Vert Tail 1 to CDo 
    end

    %% Vertical Tail #2 Contribution to CDo
    if Design_Input.Swet_v2(n)~=0 % If this component exists:
        Re_v2 = Design_Input.MAC_v2(n)*Design_Input.V_o(n)/ATMOS.nu(n); %Horz Tail Re
        Cf_v2 = .074/(Re_v2^.2); %Flat Plate Coef of Friction for Turbulent Flow
        FF_v2(n) = (1+ .6/Airfoil.X_thick_v2(n)*Airfoil.Thick_v2(n)+100*Airfoil.Thick_v2(n)^4)*...
        (1.35*M^.18*cosd(Design_Input.Sweep_v2(n))^.28); %Horz Tail Form Factor
        CDo_v2(n) = Cf_v2*FF_v2(n)*Design_Input.Q_v2(n)*Design_Input.Swet_v2(n); %Contribution of Vert Tail 2 to CDo
    end

    %% Misc. and L&P Contributions to CDo
    if Design_Input.Abase_f(n)~=0 % If this component exists:
        D_q_base = (.139+.419*(M-.161)^2)*pi*(.055)^2
        CDo_misc(n) = D_q_base/Design_Input.Sref_w(n); %Contribution of base drag to CDo
    end
    
    %% Leakage and Proturbance Contribution to CDo
    PercentLeakage = 15; % Your choice depending on how bad you think your fabrication will be 
    CDo_lp(n) = CDo_f*PercentLeakage/100; %Increase in parasite drag due to leakage and protuberance usually 3-15% of total CDo, but here just taking fuselage contribution
 
    %%Total Parasite Drag and Wetted Area
    Swet_tot(n) = Design_Input.Swet_w(n)+Design_Input.Swet_f(n)+Design_Input.Swet_h1(n)+Design_Input.Swet_h2(n)+...
        Design_Input.Swet_v1(n)+Design_Input.Swet_v2(n); %Total Wetted Area
    CDo(n) = ((CDo_w+CDo_f+CDo_v1+CDo_v2+CDo_h1+CDo_h2)/(Design_Input.Sref_w(n)))+CDo_misc+CDo_lp(n); %Total Parasite Drag Coefficient

% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
end
    
%% Oraganize into tables for output
Parasite_Data = table(CDo, CDo_f, CDo_w, CDo_h1, CDo_h2, CDo_v1, CDo_v2, CDo_misc, CDo_lp, Swet_tot);
FF_Table = table(FF_f, FF_w, FF_h1, FF_h2, FF_v1, FF_v2);

end

