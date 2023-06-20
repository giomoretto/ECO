function [kappa,specR,thetaCyl,xiO2,xBG,xBZ,zetaComb] = StaticCylinderConditions(QComb,pCyl,VCyl,parModel,parOP,NumSym)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Derive In-Cylinder Gas Mixture Species
   species = {'N2','O2','CO2','H2O','CxHy'};
 % Change in species masses and mass fractions due to combustion based
 % on global 1 step reaction
   % C_x H_y + (x + y/4)*O2 + N2 = x*CO2 + y/2*H2O + N2
   % --> nuN2  = 0
   % --> nuO2  = (-1)*(x + y/4) 
   % --> nuCO2 = x
   % --> nuH2O = y/2
   % Conversion rate of fuel at time step k in [mol/°CA]
 % Calculate Burned Fuel Mass and Mols
   MFuelConv   = QComb / parModel.thermo.fuel.lowHeatVal;
   MolFuelConv = MFuelConv / parModel.thermo.gas.molarMass.CxHy;
 % Calculate species masses [kg]
   MCylTot = 0;
 % Initialization of MCyl necessary for simulink
   for s = 1:length(species)
       MCyl.(species{s}) = 0;
   end
   MCyl.air = 0;
   MCyl.bg  = 0;
   for s = 1:length(species)-1
       MCyl.(species{s}) = ...
           parModel.thermo.fuel.nu.(species{s})*MolFuelConv*...
           parModel.thermo.gas.molarMass.(species{s}) + ...
           parOP.mCylTot * parOP.(['xi',species{s}]);
       MCylTot = MCylTot + MCyl.(species{s});
   end
   MCyl.(species{end}) = 0;
 % Calculate air and burnt gas share (should be a function of CO2)
   MCyl.air =...
       parModel.thermo.fuel.nu.air*MolFuelConv*...
       parModel.thermo.gas.molarMass.air + ...
       parOP.mCylTot * parOP.xiAir;
   MCyl.bg = MCylTot - MCyl.air;
 % Calculate Species Mass Fractions [-]
   for s = 1:length(species)
       xi.(species{s}) = MCyl.(species{s}) ./ MCylTot; 
   end
   xi_air = MCyl.air / MCylTot;
   xi_bg  = 1 - xi_air;
 % Calculate Specific Gas Constant [J/(kg*K)]
   specR = 0;
   for k = 1:length(species)
       specR = specR + xi.(species{k})*...
           parModel.thermo.gas.specR.(species{k});
   end
    
%% Derive Temperature and Kappa
 % Temperature [°K]
   thetaCyl = pCyl.*VCyl./(MCylTot.*specR);
 % Kappa update
   if strcmp(NumSym,'Sym')
       kappa = calcKappa_casADi(xi,thetaCyl,parModel.thermo.gas.specR);
   else
       kappa = calcKappa(xi,thetaCyl,parModel.thermo.gas.specR);
   end
 % oxygen mass fraction   
  xiO2 = xi.(species{2});
  xBG  = xi_bg;

 % burned zone mass fraction
 xBZ = QComb*(parModel.thermo.fuel.airFuelRatioSt/(1-parOP.xiBg)+1)/...
     (MCylTot*parModel.thermo.fuel.lowHeatVal+QComb);

end

