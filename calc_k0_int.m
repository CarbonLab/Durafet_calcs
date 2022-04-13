function [k0, Eoinsitu] = calc_k0_int(Vint, TC, pH_insitu, UseDefaultk2)

% [Eoinsitu, Eo20C] = calc_Eo_int(Vint, TC, pH_insitu)
%
% Calculate Eo (calibration coefficients) for Durafet using the internal
% reference electrode at in situ and 25C. 
%
% Created by: Yui Takeshita
% Monterey Bay Aquarium Research Institute
% Version 1 Created: November 23, 2016


R = 8.31451; %Universal Gas Constant
F = 96487; %Faraday Constant

if(strcmp(UseDefaultk2, 'UseDefaultk2'))
    k2 = -0.001455; %From Martz et al. 2010
else
    k2 = UseDefaultk2;
end

    
%Calcuate Eo
TK = TC + 273.15;
S = R.*(TK).*log(10)/F;
Eoinsitu = Vint - S.*pH_insitu;

k0 = Eoinsitu - k2.*TC;

return
