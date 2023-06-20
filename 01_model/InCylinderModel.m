function [xdot,y] = InCylinderModel(x,u,ca,parModel,parOP,NumSym)
%In-Cylinder Model which calculates the state derivatives, given the states
%and the input.
%   States x       = [p [Pa], QComb [J], IMEP [Pa]]
%   input  u       = dQcomb [J/degCA]
%   output y       = [thetaCyl [°K] kappa [-] dQWall [J/degCA] specR [J/(kg*K)]]
%   Crank Angle    = ca [degCA]
%   Disturbance d  = w_e (engine speed, [1/s])
%   State Derivatives xdot are gradients with respect to degCA [stateUnit/degCA]


%% Define States
    pCyl  = x(1);
    QComb = x(2);
    % IMEP  = x(3);
  
%% Derive Actual Volume and Volume Derivative
 % VCyl: [m^3] / dVCyl: [m^2/degCA] / stroke: [m]
   [VCyl,dVCyl,stroke] = CylVol(ca,parModel.eng);
  
%% Derive kappa, specific gas constant and Temperature
   [kappa,specR,thetaCyl] = StaticCylinderConditions(QComb,pCyl,VCyl,...
                                                    parModel,parOP,NumSym);


%% Derive Heat Release rate [J/degCA]
   dQComb = u;

%% Derive In-Cylinder Derivatives
% % Wall Heat Transfer [J/degCA], Cylinder Pressure Derivative [Pa/degCA]
% % and Indicated Mean Effective Pressure Derivative [Pa/degCA]
    [dpCyl,dIMEP,dQWall] = DynamicCylinderConditions(ca,pCyl,dQComb,...
        VCyl,dVCyl,stroke,thetaCyl,kappa,specR,parModel,parOP);
   
%% Write Derivatives to xdot
   xdot = [dpCyl; dQComb; dIMEP];
   y    = [thetaCyl;kappa;dQWall;specR];

end
