function [pH_int] = calc_dfet_pHint_k0(Vint, tempC, k0, UseDefaultk2)

% pH_int = calc_dfet_pHint(Vint, Vtherm, Vint_calib, Vtherm_calib, pH_calib)
%
% calculates pH from the durafet using the internal reference electrode. 
%
% Created by: Yui Takeshita
% Monterey Bay Aquarium Research Institute
% Version 1 Created: November 23, 2016

R = 8.31451; %Universal Gas Constant
F = 96487; %Faraday Constant

% decide whether to use 
if(strcmp(UseDefaultk2, 'UseDefaultk2'))
    k2 = -0.001455; %From Martz et al. 2010
else
    k2 = UseDefaultk2;
end


tempK = tempC + 273.15;
%Calculate pH accoriding to Nernst Equation
S = R*tempK*log(10)/F; %Nernst Slope at given temperature
Eoint = k0 + k2.*tempC; %Eo at given temperature
pH_int = (Vint - Eoint)./S;

