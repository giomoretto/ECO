function [tau_chem] = ignitionDelayJoerg(p,T,XO2)
%Ignition delay identified by Joerg
% - p: cylinder pressure in [Pa]
% - T: cylinder temperature in [K]
% - XO2: mass fraction oxygen in [0-1] 
%   Detailed explanation goes here
pRef = 40e5;
XO2ref = 0.224;
R = 1;

fuelType = 'Diesel';

funcArr = @(p,T,XO2,cPres,cO2,k,Ea) k.*(p./pRef).^(cPres).*(XO2./XO2ref).^(cO2).*exp(Ea./T./R);


switch fuelType
    case 'Diesel'
        tau_HT = funcArr(p,T,XO2,-1.05,-1.5,1.25e-07,15813);
        tau_LT = funcArr(p,T,XO2,0    ,0   ,1.09e-08,13188);
        
        A_HT   = funcArr(p,T,XO2,-1   ,-1  ,3.11e-03,0);
        A_LT   = funcArr(p,T,XO2,0    ,0   ,4.51e-12,21079);
        
    case 'nHeptane'
        tau_HT = funcArr(p,T,XO2,-1.1,-1.5,1.3e-07,15813);
        tau_LT = funcArr(p,T,XO2,0   ,-0.2,2.65e-08,13188);
        
        A_HT   = funcArr(p,T,XO2,-0.5,-0.7,3.39e-3,0);
        A_LT   = funcArr(p,T,XO2,0    ,1,4.51e-12,21079);
        
end

w_HT = A_HT./(A_HT+A_LT);
w_LT = 1-w_HT;

tau_chem = w_HT.*tau_HT + w_LT.*tau_LT;


end

