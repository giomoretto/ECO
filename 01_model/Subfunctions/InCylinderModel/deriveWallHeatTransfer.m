function [dQWalldPhi] = deriveWallHeatTransfer(ca,stroke,...
                                             thetaCyl,pCyl,...
                                             thetaIvc,pIvc,volIvc,vol,...
                                             engSpd,parModel)

% Transformation factor dPhidT in [degCA/s]
dPhidT = engSpd*360;

% Transformation factor dTdPhi in [s/degCA]
dTdPhi = 1/dPhidT;

% Wall temperature in [K]
airExcessFactor = 1;
thetaWall = 360 + 9*airExcessFactor^0.4*sqrt(parModel.eng.bore*engSpd*60) ;%+ Par.wall.thetaOffset;
% thetaWall = 440;

% Surface area of the combustion chamber, area  of top/fire land not
% considered and area of piston assumed to be equal to bore area
areaCombChamber =...
    parModel.eng.areaCylHead + parModel.eng.areaBore + pi*parModel.eng.bore*stroke;

% Mean piston velocity cm in [m / s]
cm = 4*parModel.eng.radCrank*engSpd;

% Constant c1
c1 = parModel.wall.c1;

% Characteristic speed c in [m / s]
% if phi >= soi
    % Estimated motored pressure
    pCylMot = pIvc*(volIvc/vol)^1.33;
    
    deltaPCyl = pCyl - pCylMot;
    
%     if thetaWall < 600
%         c2 = 3.24e-3;
%     else
        c2 = 5.0e-3 + 2.3e-5*(thetaWall - 600);
%     end
    
    c = c1*cm + c2*(parModel.eng.volDis*thetaIvc)/(pIvc*volIvc)*deltaPCyl*parModel.wall.scaleFiredWHL;
%     c = c1*cm;

% else
%     c = c1*cm;
% end

% Heat transfer coefficient in [J / s*m^2*K]
alpha = 0.013*parModel.eng.bore^(-0.2)*pCyl^(0.8)*thetaCyl^(-0.53)*c^(0.8);
% scaling identified with compression curve
alpha = alpha*parModel.wall.scaleWHL;
% enable ca variable wall heat losses to obtain correct center of comb.
if parModel.Opts.EnableVariableWHL
    
    scaleWHL_CAvariable = caVariable_WHL(ca,parModel);
    alpha = alpha*scaleWHL_CAvariable;
    
end



% Heat flow through cylinder wall in [J/s]
dQWalldT = alpha*areaCombChamber*(thetaWall - thetaCyl);

% Heat flow through cylinder wall in [J/degCA]
% dQWalldPhi = 0.5*dQWalldT*dTdPhi;
dQWalldPhi = 0.8*dQWalldT*dTdPhi;


end

