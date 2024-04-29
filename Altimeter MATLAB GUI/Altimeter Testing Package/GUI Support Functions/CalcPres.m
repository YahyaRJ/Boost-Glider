function [Pressure, Temp] = CalcPres(C1,C2,C3,C4,C5,C6,P,T)
%CALCPRES Calculates pressure and temperature.
%   [Pressure, Temp] = CalcPres(C1,C2,C3,C4,C5,C6,P,T) uses calibration
%   coefficients C1, C2, ... , C6 to convert raw pressure and temperature
%   values from the ASEN 2004 Altimeter Rev.9 using pressure sensor MS5611
%   to real pressure and temperature values in mBarr and degrees C.
%
%   See data sheet for MS5611 Pressure Sensor pg.7
%http://www.meas-spec.com/downloads/MS5611-01BA03.pdf (July 2014)
dT = T - C5*2^8;
Temp = 2000+dT*C6/2^23;
Off = C2*2^16+(C4*dT)/2^7;
Sens = C1*2^15+(C3*dT)/2^8;
Pressure = ((P*Sens/2^21-Off)/2^15)/100;
Temp = Temp/100;
end
