function [Par] = ModelParameters(varargin)

%% Engine geometry

  % Total displacement volume in [m^3]
    eng.volDisTotal = 1.968e-3;
  % Number of cylinders in [1]
    eng.numCyl = 4;
  % Displacement volume per cylinder in [m^3]
    eng.volDis = eng.volDisTotal/eng.numCyl;
  % Compression ratio (epsilon) in [1]
    eng.epsilon = 16.5;
  % Clearance volume per cylinder in [m^3]
    eng.volClear = eng.volDis/(eng.epsilon - 1);
  % Stroke in [m]
    eng.stroke = 95.5e-3;
  % Crank radius in [m]
    eng.radCrank = eng.stroke/2;
  % Con-rod length in [m]
    eng.lengthConRod = 3*eng.radCrank;
  % Bore area in [m^2]
    eng.areaBore = eng.volDis/eng.stroke;
  % Cylinder head area in [m^2]
    eng.areaCylHead = eng.areaBore*1.2;
  % Bore diameter in [m]
    eng.bore = 2*sqrt(eng.areaBore/pi);

%% Thermodynamic fuel and gas properties

% % Species amount fractions for air in [1]
   thermo.psi.air.O2 = 0.21;
   thermo.psi.air.N2 = 0.79;

% % Fuel related parameters (Diesel)
  % Lower heating value in [J/kg] (Source: Merker et al. 2011)
    thermo.fuel.lowHeatVal = 42.6e6;
  % Upper heating value in [J/kg]
    thermo.fuel.upHeatVal = 45e6;
  % C and H composition, C_x H_y
    thermo.fuel.x = 10.8; 
    thermo.fuel.y = 18.7;
  % Stochiometric coefficients for global combustion reaction of fuel surrogate
  % with oxygen assuming complete combustion  
  % CxHy + (x + y/4)*O2 = x*CO2 + y/2*H20 
    thermo.fuel.nu.N2   = 0;
    thermo.fuel.nu.O2   = (-1)*(thermo.fuel.x + thermo.fuel.y/4);
    thermo.fuel.nu.CO2  = thermo.fuel.x;
    thermo.fuel.nu.H2O  = thermo.fuel.y/2;
    thermo.fuel.nu.CxHy = -1;
    % CxHy + 1/xiAirO2*(x + y/4)*air = 1*bg
    thermo.fuel.nu.air = (-1)*1/thermo.psi.air.O2*...
                         (thermo.fuel.x + thermo.fuel.y/4);
    thermo.fuel.nu.bg  = 1;
  % Density at 22°C in [kg/m^3] 
    thermo.fuel.rho = 834*(1+0.00067*7)*1;
  % Kinematic fuel viscosity in [m^2/s]
    thermo.fuel.kinVisc = 2.8e-6;
  % Surface tension in [N/m]
    thermo.fuel.surfTen = 0.03;
  % Evaporation temperature in [K]
    thermo.fuel.thetaEvap = 570;
  % Evaporation enthalpy in [J/kg]
    thermo.fuel.evapEnthalpy = 359*1e3;
  % Spec. heat capacity in [J/(kg K)]
    thermo.fuel.specHeatCap = 1926;
    
% % In-Cylinder Gas related parameters
  % Species amount fractions for burnt gas in [1]
    molBgTot = thermo.fuel.x + thermo.fuel.y/2 +...
        thermo.psi.air.N2/thermo.psi.air.O2*(thermo.fuel.x + thermo.fuel.y/4);
    thermo.psi.bg.H2O = (thermo.fuel.y/2)/molBgTot;
    thermo.psi.bg.CO2 = thermo.fuel.x/molBgTot;
    thermo.psi.bg.N2  = (thermo.psi.air.N2/thermo.psi.air.O2*...
        (thermo.fuel.x + thermo.fuel.y/4))/ molBgTot;
  % Species Molar masses in [kg/mol]
    thermo.gas.molarMass.C    = 12.0108e-3;
    thermo.gas.molarMass.H    = 1.0079e-3;
    thermo.gas.molarMass.O    = 15.9994e-3;
    thermo.gas.molarMass.N    = 14.0067e-3;
    thermo.gas.molarMass.O2   = 2*thermo.gas.molarMass.O;
    thermo.gas.molarMass.N2   = 2*thermo.gas.molarMass.N;
    thermo.gas.molarMass.H2O  = 2*thermo.gas.molarMass.H + ...
                                thermo.gas.molarMass.O;
    thermo.gas.molarMass.CO2  = thermo.gas.molarMass.C + ...
                                2*thermo.gas.molarMass.O;
    thermo.gas.molarMass.air  = thermo.psi.air.O2*thermo.gas.molarMass.O2+ ...
                                thermo.psi.air.N2*thermo.gas.molarMass.N2;
    thermo.gas.molarMass.fuel = thermo.fuel.x*thermo.gas.molarMass.C + ...
                                thermo.fuel.y*thermo.gas.molarMass.H;
    thermo.gas.molarMass.CxHy = thermo.fuel.x*thermo.gas.molarMass.C + ...
                                thermo.fuel.y*thermo.gas.molarMass.H;
    thermo.gas.molarMass.bg   = thermo.psi.bg.H2O*thermo.gas.molarMass.H2O +...
                                thermo.psi.bg.CO2*thermo.gas.molarMass.CO2 +...
                                thermo.psi.bg.N2*thermo.gas.molarMass.N2;
  % Universal/molar gas constant in [J/(K*mol)] 
    thermo.gas.univR = 8.314;
  % Specific gas constants in [J/K*kg]
    thermo.gas.specR.O2   = thermo.gas.univR/thermo.gas.molarMass.O2;
    thermo.gas.specR.N2   = thermo.gas.univR/thermo.gas.molarMass.N2;
    thermo.gas.specR.H2O  = thermo.gas.univR/thermo.gas.molarMass.H2O;
    thermo.gas.specR.CO2  = thermo.gas.univR/thermo.gas.molarMass.CO2;
    thermo.gas.specR.air  = thermo.gas.univR/thermo.gas.molarMass.air;
    thermo.gas.specR.fuel = thermo.gas.univR/thermo.gas.molarMass.fuel;
    thermo.gas.specR.CxHy = thermo.gas.univR/thermo.gas.molarMass.CxHy;
    thermo.gas.specR.bg   = thermo.gas.univR/thermo.gas.molarMass.bg;
  % Species mass fractions for air in [1]
    thermo.xi.air.O2 = thermo.psi.air.O2*thermo.gas.molarMass.O2/...
                       thermo.gas.molarMass.air;
    thermo.xi.air.N2 = thermo.psi.air.N2*thermo.gas.molarMass.N2/...
                       thermo.gas.molarMass.air;
  % Species mass fractions for burnt gas [1]
    thermo.xi.bg.CO2 = thermo.psi.bg.CO2*thermo.gas.molarMass.CO2/...
                       thermo.gas.molarMass.bg;
    thermo.xi.bg.H2O = thermo.psi.bg.H2O*thermo.gas.molarMass.H2O/...
                       thermo.gas.molarMass.bg;
    thermo.xi.bg.N2  = thermo.psi.bg.N2*thermo.gas.molarMass.N2/...
                       thermo.gas.molarMass.bg;
  % Stoichiometric air fuel ratio (air requirement) [kg air/kg fuel] = [1]
    thermo.fuel.airFuelRatioSt = (1/thermo.psi.air.O2*...
        (thermo.fuel.x + thermo.fuel.y/4)*thermo.gas.molarMass.air)/...
        (thermo.gas.molarMass.CxHy);

%% Wall heat loss model (Woschni)
% enable wall heat loss model
  Opts.EnableWallHeatLoss = true;
% enable variable wall heat loss
  Opts.EnableVariableWHL  = false;
% Intake swirl number in [1]
  wall.intSwirlNumber = 2.5;
% C1 factor in [1]
  wall.c1 = 2.28 + 0.308*wall.intSwirlNumber;
% scale factor
  wall.scaleWHL = 2.8;
  wall.scaleFiredWHL = 0.05;
% ca-variable wall heat losses
  wall.smoothness = 10;
  wall.dropPercentage = 25;
  wall.caDropLocation = 0;



%% Ignition Delay Model (Simple)

% % select ignition delay model
    Opts.IgnitionDelayModel = 'Joerg'; % 'Joerg', 'Nils'

% % Physical ignition delay modeled as constant in [s]
    ign.phys.const = 2.5e-4;

% % minimal and maximal ignition delay
    ign.tauMin = 0; %0
    ign.tauMax = 2e-3; %2e-3

% % correction factor for physical ignition delay
    ign.chem.const = 1;
%% Combustion heat release rate model parameters (simple)

% Parameters for Simplified Model
    comb.cDif.a           = 5.27e-4;
    comb.cDif.b           = 2.07e3;
    comb.cPre.a           = 0;
    comb.cPre.b           = 3.8e3;
    comb.gamma.tauMin     = 2.7e-4;
    comb.gamma.tauMax     = 2.9e-4; 
% Combustion efficiency in [1]
    comb.etaComb          = 1;



%% Injector model parameters
    inj = struct(); % write algebraic model parameters here
    

%% Gather all parameters in overall structure
    Par        = struct();
    Par.eng    = eng;
    Par.thermo = thermo;
    Par.wall   = wall;
    Par.ign    = ign;
    Par.inj    = inj;
    Par.comb   = comb;
    Par.Opts   = Opts;
    
    
    if nargin == 1
        xSol = varargin{1};
        Par = loadParameters(Par,xSol);
    end

end


function [parModel] = loadParameters(parModel,xSol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Simulation With optimal Combustion
parModel.comb.cDif.a = xSol(1);
parModel.comb.cDif.b = xSol(2);
parModel.comb.cPre.a = xSol(3);
parModel.comb.cPre.b = xSol(4);
parModel.ign.phys.const = xSol(5);
parModel.comb.gamma.tauMin = xSol(6);
parModel.comb.gamma.tauMax = xSol(7);
parModel.comb.ign.tauMin = xSol(8);
parModel.comb.ign.tauMax = xSol(9);

end



