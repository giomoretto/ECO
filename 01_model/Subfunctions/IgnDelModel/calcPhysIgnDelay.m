function [tauPhys] = calcPhysIgnDelay(thetaCyl,pCyl,specR,kappa,pRail,Par,varargin)

% Default value
physIgnDelayModel = 1;

if nargin > 6
    physIgnDelayModel = varargin{1};
else
end

switch physIgnDelayModel
    case 1      % Constant value
        tauPhys = Par.ign.phys.const;
        
    case 2      % Joerg 2018, Diss., VKA, RWTH Aachen University
        % Estimation for pressure drop over nozzle in [Pa]
        deltaPNoz = pRail - pCyl;

        % Density of in-cylinder charge in [kg/m^3]
        rhoCyl = pCyl/(specR*thetaCyl);

        % Fuel jet break-up time in [s]
        tauBreakUp =...
            4.351*Par.ign.phys.dNozzle*Par.thermo.fuel.rho/...
            (Par.ign.phys.disCoef^2*sqrt(rhoCyl*deltaPNoz));

        % Sauter mean diameter (droplet volume to surface ratio) in [m]
        smd =...
            3.08*...
            Par.thermo.fuel.kinVisc^(0.385)*...
            Par.thermo.fuel.surfTen^(0.737)*...
            Par.thermo.fuel.rho^(0.737)*...
            rhoCyl^(0.06)*...
            deltaPNoz^(-0.54);

        % Isobaric specific heat capacity of the in-cylinder charge in [J/kgK]
        specHeatCap = specR*(kappa/(kappa - 1));

        % Fuel evaporation rate in [m^2/s]
        fuelEvapRate =...
            8*Par.ign.phys.gas.thermCond/(rhoCyl*specHeatCap)*...
            log(1 + specHeatCap/Par.thermo.fuel.evapEnthalpy*(thetaCyl - Par.thermo.fuel.thetaEvap));

        % fuelEvapRate =...
        %     8*Par.ign.phys.gas.thermCond/(rhoCyl*Par.thermo.fuel.specHeatCap)*...
        %     log(1 + Par.thermo.fuel.specHeatCap/Par.thermo.fuel.evapEnthalpy*(thetaCyl - Par.thermo.fuel.thetaEvap));

        % Evaporation time in [s]
        tauEvap = 1/fuelEvapRate*smd^2;

        % Total physical ignition delay time in [s]
        tauPhys = tauBreakUp + tauEvap;
        
    otherwise
        tauPhys = 0;
        
end



end

