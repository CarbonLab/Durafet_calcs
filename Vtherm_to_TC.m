
function [tempC,Rtherm] = Vtherm_to_TC(Vtherm)

%This function takes the Vtherm from the durafet and calculates temperature
%in celcius. Input can be a vector. 

T1 = 340.98198629984;
T2 = -0.0000910256997011314;
T3 = -95.088066665817;   
T4 = 0.965370273924651;
Rtherm =20000./(3.3./Vtherm-1);
tempC =  T1 + T2*Rtherm + T3*log10(Rtherm) + T4*(log10(Rtherm)).^3;

return