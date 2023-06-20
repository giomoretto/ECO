function [dpCyl,dIMEP,dQWall] = DynamicCylinderConditions(ca,pCyl,dQComb,...
    VCyl,dVCyl,stroke,thetaCyl,kappa,specR,parModel,parOP)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% % Derive Wall Heat Transfer [J/degCA]
    if parModel.Opts.EnableWallHeatLoss
        dQWall = deriveWallHeatTransfer(ca,stroke,thetaCyl,pCyl,...
                                        parOP.thetaIvc,parOP.pInt,...
                                        parOP.VInt,VCyl,parOP.engSpd,parModel);
    else
        dQWall = 0;
    end

% % Derive Cylinder Pressure Derivative [Pa/degCA]
  % Get fuel massflow derivative assuming instantaneous combustion
%     dMFuel   = 0; % dQComb / Par.thermo.fuel.lowHeatVal;
    dMFuel = 0;
    dpCyl = -kappa * pCyl * dVCyl / VCyl + ...
            specR*thetaCyl*dMFuel / VCyl + ...
            (kappa - 1) / VCyl * (dQWall + dQComb);
        
   
% % Derive Indicated Mean Effective Pressure Derivative [Pa/degCA]
    dIMEP = pCyl*dVCyl/parModel.eng.volDis;


end

