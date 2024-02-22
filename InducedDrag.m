function [InducedDrag_Data] =...
    InducedDrag(Design_Input,WingLiftModel,WingLiftCurve,WingDragCurve,WingGeo_Data,Count,Airfoil,Parasite_Data)
%% Induced Drag Model Function Summary
% This function evaluates different Oswalds Efficiency Factor models for 
% use in your drag polar model.  It compiles and outputs a variables table
% for three different Oswalds models (mod1, mod2, mod3).  For each Oswalds
% model, the Oswalds Efficiencty Factor (eo_mod1,2,3), and the k1 values
% (k1_mod1,2,3) are outputted in the InducedDrag_Data table.  Note that
% depending on the Oswalds models chosen to evaluate, you may or may not
% need information from teh WingGeo_Data table from the WingGeo function.
% Additionally, this code supports the calculation of the k2 values for
% evaluating non-symmetric airfoil design, but is not required.

%% Outputs:
%
% InducedDrag_Data:
%   Table containing lift curve slope, zero lift AoA, and oswalds info
%   (columns) for each input (rows)


%% Preallocate variables of interest
eo_mod1 = zeros(Count, 1); % Oswalds for Model #1
eo_mod2 = zeros(Count, 1); % Oswalds for Model #2
eo_mod3 = zeros(Count, 1); % Oswalds for Model #1
k1_mod1 = zeros(Count, 1); % k1 for Model #1
k1_mod2 = zeros(Count, 1); % k1 for Model #2
k1_mod3 = zeros(Count, 1); % k1 for Model #3
k2_mod1 = zeros(Count, 1); % k2 for Model #1
k2_mod2 = zeros(Count, 1); % k2 for Model #2
k2_mod3 = zeros(Count, 1); % k2 for Model #3

% NOTE: k2 values not required if only symmetric airfoils used; however,
% this version of the code includes it as an option

%% Loop through different configurations
for n = 1:Count 
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
    %Find CL min Drag value of wing drag polar to estimate k2
    [CD_min,CD_min_index] = min(WingDragCurve{n,:});
    CL_minD = WingLiftCurve{n,CD_min_index};

    %Cavallo Oswalds Model (Baseline required)
    eo_mod1(n) = 1.78*(1-0.045*Design_Input.AR_w(n)^0.68)-0.64;
    k1_mod1(n) =  1/(pi*eo_mod1(n)*Design_Input.AR_w(n));
    k2_mod1(n) = -2*k1_mod1(n)*CL_minD;

    %Student Option 1 Oswalds Model (Obert)

    eo_mod2(n) = 1/(1.05+0.007*pi*Design_Input.AR_w(n)); %Oswalds Estimate
    k1_mod2(n) = 1/(pi*eo_mod2(n)*Design_Input.AR_w(n)); 

   
    %Student Option 2 Oswalds Model (McCormick)

    eo_mod3(n) = 1/(1+0.03); %Oswalds Estimate
    k1_mod3(n) = 1/(pi*eo_mod3(n)*Design_Input.AR_w(n));

 
    %Student Option 3 Oswalds Model (Grosu)

    eo_mod4(n) = 1/(1.08+(0.028*Airfoil.Thick_w(n))/(CL_minD^2)*pi*Design_Input.AR_w(n)); %Oswalds Estimate
    k1_mod4(n) = 1/(pi*eo_mod1(n)*Design_Input.AR_w(n));

    %Student Option 4 Oswalds Model (Kroo)
    s = 1-2*(Design_Input.Dia_f(n)/WingGeo_Data.b_w(n))^2;
    Q = 1/(0.99*s);
    P = 0.38*Parasite_Data.CDo(n); 
    eo_mod5(n) = 1/(Q+P*pi*Design_Input.AR_w(n)); %Oswalds Estimate
    k1_mod5(n) = 1/(pi*eo_mod1(n)*Design_Input.AR_w(n));

     %Student Option 5 Oswalds Model (Schaufele)

    eo_mod6(n) = 1/(1.03+0.38*Parasite_Data.CDo(n)); %Oswalds Estimate
    k1_mod6(n) = 1/(pi*eo_mod1(n)*Design_Input.AR_w(n)); 

% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////   
end

%% Oraganize into table for output
InducedDrag_Data = table(eo_mod1, eo_mod2, eo_mod3, k1_mod1, k1_mod2, k1_mod3, k2_mod1, k2_mod2, k2_mod3);

end


