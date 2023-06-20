function [xdot,y] = CompleteModel(x,u,ca,parModel,parOP,NumSym)
%In-Cylinder Model which calculates the state derivatives, given the states
%and the input.
%   States x       = [1)pCyl [Pa], 2)QComb [J], 3)IMEP [Pa], 4)thetaUZ [K], 5)NO [mole fraction in BG]
%   input  u       = [SOE [degCA ATDC] DOE [mu s]]
%   output y       = [1) thetaCyl [°K] 2)kappa [-] 3)dQWall [J/degCA] 4)specR [J/(kg*K)] ...
%                     5)mFuelPrep [] 6)tauIgn [] 7)xBG [] 8)thetaBZ [K] 9)NOppm [ppm whole gas] ...
%                     20)Phi [] 11)NOeq[mol/cm3] 12)xBZ []] 
%   Crank Angle    = ca [degCA]
%   Disturbance d  = w_e (engine speed, [1/s])
%   State Derivatives xdot are gradients with respect to degCA [stateUnit/degCA]

% % Sanity checks
  % Limit states to physical boundaries
    x   = checkValidity(x,ca,parModel,parOP);

% % Define Inputs and States
  % Inputs
    nInj    = numel(u)/2;
    SOE     = u(1:nInj);
    DOE     = u(nInj+1:2*nInj);
  % States
    pCyl    = x(1);
    QComb   = x(2);
    % IMEP  = x(3);
    if numel(x) > 3
        thetaUZ = x(4);
        NO      = x(5);
    end

% % Derive Actual Volume and Volume Derivative
  % VCyl: [m^3] / dVCyl: [m^3/degCA] / stroke: [m]
    [VCyl,dVCyl,stroke] = CylVol(ca,parModel.eng);
  
% % Derive kappa, specific gas constant and Temperature
    [kappa,specR,thetaCyl,xiO2,xBG,xBZ] = StaticCylinderConditions(QComb,pCyl,...
        VCyl,parModel,parOP,NumSym);
    
% % Derive Ignition Delay [s]
   tauIgn = ignDelModel(pCyl,thetaCyl,xiO2,parModel);

% % Derive Prepared Fuel Mass [kg]
  % Actual Crank angle minus ignition delay, converted to microseconds
  % after top dead centre
    tPrep  = (ca/(parOP.engSpd*360)-tauIgn)*1e6;
  % Start of Energyzing converted to microseconds after to dead centre
    tSOE  = SOE./(parOP.engSpd*360)*1e6;
  % Injected fuel at actual crank angle minus ignition delay [kg]
    [mFuelPrep,mTot,~] = AlgebraicInjectorModel(tPrep,parOP.pRail,tSOE,DOE);
    
  % calculate Phi = 1/lambda. Phi is used in opt. to prevent unrealisitic
  % amounts of injected fuel
    Phi = mFuelPrep*parModel.thermo.fuel.airFuelRatioSt*parModel.thermo.xi.air.O2...
        / parOP.mCylTot/parOP.xiO2;

% % Derive Heat Release Rate [J/degCA]
    dQComb = CombustionModel(QComb,mFuelPrep,tauIgn,parOP.engSpd,parModel);
      
% % Wall Heat Transfer [J/degCA], Cylinder Pressure Derivative [Pa/degCA]
% % and Indicated Mean Effective Pressure Derivative [Pa/degCA]
    [dpCyl,dIMEP,dQWall] = DynamicCylinderConditions(ca,pCyl,dQComb,...
        VCyl,dVCyl,stroke,thetaCyl,kappa,specR,parModel,parOP);

if numel(x) > 3
% % Derive dThetaUZ and thetaBZ [K] from a Two-Zone Model
    [dThetaUZ,thetaBZ] = TwoZoneModel(QComb,thetaUZ,dpCyl,pCyl,dVCyl,VCyl,...
        thetaCyl,kappa,mTot,parOP,parModel,xBG,xBZ);
    
% % Derive NOx concentration in the exhaust
    [dNO,NOppm,NOeq] = NOxModel(QComb,thetaUZ,thetaBZ,NO,mTot,parOP,parModel,...
        xBG,VCyl,dVCyl,pCyl,dpCyl,thetaCyl);
end

% % Write Derivatives to xdot
    if numel(x) > 3
        xdot = [dpCyl; dQComb; dIMEP; dThetaUZ; dNO];
        y    = [thetaCyl;kappa;dQWall;specR;mFuelPrep;tauIgn;xBG;thetaBZ;NOppm;Phi;NOeq;xBZ];
    else
        xdot = [dpCyl; dQComb; dIMEP];
        y    = [thetaCyl;kappa;dQWall;specR;mFuelPrep;tauIgn;Phi;xiO2];
    end
    
    
end
