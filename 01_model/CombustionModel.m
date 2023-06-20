function [dQCombdPhi] = CombustionModel(qComb,mFuelPrep,tauIgn,we,Par)
%COMBUSTIONMODEL Is taken from Nils master's thesis and is based  on a
%Chmela combustion model, mostly taken from Jörg 2018 (Diss)
%   Inputs: qComb:      total heat released by burned fuel so far [J]
%           mFuelPrep:  total injected fuel, which is ready to be burned,
%                       hence mFuel(t-tauIgn) [kg]
%           tauIgn:     Ignition Delay [s]
%           we:         Engine Speed [1/s]
%           Par:        Parameter Set
%  Outputs: dQCombdPhi: Heat Release Rate in [J/degCA]

    % time passed per °CA and vice versa
    dPhidT = we*360;
    dTdPhi = 1/dPhidT;

    % Available fuel energy potential in [J]
    qDis      = 0; % General Assumption in Nils Models
    dQDisdPhi = 0;
    qFuelPot = (mFuelPrep*Par.thermo.fuel.lowHeatVal*Par.comb.etaComb - (qComb - qDis));

    % Limit qFuelPot to 0
%     eps = 1e-3; qFuelPot = 1/2*qFuelPot + 1/2*sqrt(eps^2 + qFuelPot^2);
    
    % Simplified version of empirical mixing rate modeling
    % approach (C = a*dQ/dt + b)
    % Weighting factor for mixing rate in [1]
    % gamma = 0    : pure diffusive mixing rate
    % gamma = 1    : pure premixed mixing rate
    % 0 < gamma < 1: mixed mixing rate
    gamma = evalCombWeightingFun(tauIgn,...
                             Par.comb.gamma.tauMin,...
                             Par.comb.gamma.tauMax);

    % Mixing rate model parameters
    aDif = Par.comb.cDif.a;
    bDif = Par.comb.cDif.b;
    aPre = Par.comb.cPre.a;
    bPre = Par.comb.cPre.b;

    % Heat release rate in [J/degCA]
    deNom     = (1 - qFuelPot.*(gamma.*aPre + (1 - gamma).*aDif));
%     eps = 1e-5; deNom = 1/2*deNom + 1/2*sqrt(eps + deNom^2);
    
    dQCombdPhi =...
        (qFuelPot.*dTdPhi.*(gamma.*bPre +...
         (1 - gamma).*bDif) + dQDisdPhi)./...
        deNom;
    
    % Limit dQCombdPhi to 0
    eps        = 1e-4;
    dQCombdPhi = 1/2*dQCombdPhi + 1/2*sqrt(eps^2 + dQCombdPhi^2);

    % % Diffusive mixing rate in [1/s]
    % mixRate_dif = aDif*dQCombdPhi*dPhidT + bDif;

    % % Premixed mixing rate in [1/s]
    % mixRate_pre = aPre*dQCombdPhi*dPhidT + bPre;

    % % Total mixing rate in [1/s]
    % mixRate_tot = gamma*mixRate_pre + (1 - gamma)*mixRate_dif;
                
end

