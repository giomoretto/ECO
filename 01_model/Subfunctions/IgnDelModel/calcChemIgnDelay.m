function [tauChem] = calcChemIgnDelay(theta,p,xiO2,Par,varargin)

% Default value
chemIgnDelayModel = 1;

if nargin > 4
    chemIgnDelayModel = varargin{1};
else
end

switch chemIgnDelayModel
    case 1      % Single Arrhenius approach
        k1 = Par.ign.chem.k1;
        k2 = Par.ign.chem.k2;
        k3 = Par.ign.chem.k3;
        k4 = Par.ign.chem.k4;

        tauChem = k1*exp(-k2*p) + k3*exp(-k4*theta);
        
    case 2      % Joerg 2018, Diss., VKA, RWTH Aachen University
        % Reference states
        ArrPar.pRef    = Par.ign.chem.pRef;
        ArrPar.xiO2Ref = Par.ign.chem.xiO2Ref;

        % High temperature ignition delay
        ArrPar.thetaEa = Par.ign.chem.tau.highTemp.thetaEa;
        ArrPar.k       = Par.ign.chem.tau.highTemp.k;
        ArrPar.nP      = Par.ign.chem.tau.highTemp.nP;
        ArrPar.nO2     = Par.ign.chem.tau.highTemp.nO2;
        tauHigh = evalArrheniusEqn(theta,p,xiO2,ArrPar);

        % Low temperature ignition delay
        ArrPar.thetaEa = Par.ign.chem.tau.lowTemp.thetaEa;
        ArrPar.k       = Par.ign.chem.tau.lowTemp.k;
        ArrPar.nP      = Par.ign.chem.tau.lowTemp.nP;
        ArrPar.nO2     = Par.ign.chem.tau.lowTemp.nO2;
        tauLow = evalArrheniusEqn(theta,p,xiO2,ArrPar);

        % High temperature activation term
        ArrPar.thetaEa = Par.ign.chem.act.highTemp.thetaEa;
        ArrPar.k       = Par.ign.chem.act.highTemp.k;
        ArrPar.nP      = Par.ign.chem.act.highTemp.nP;
        ArrPar.nO2     = Par.ign.chem.act.highTemp.nO2;
        actHigh = evalArrheniusEqn(theta,p,xiO2,ArrPar);

        % Low temperature activation term
        ArrPar.thetaEa = Par.ign.chem.act.lowTemp.thetaEa;
        ArrPar.k       = Par.ign.chem.act.lowTemp.k;
        ArrPar.nP      = Par.ign.chem.act.lowTemp.nP;
        ArrPar.nO2     = Par.ign.chem.act.lowTemp.nO2;
        actLow = evalArrheniusEqn(theta,p,xiO2,ArrPar);

        % High and low temperature weighting terms
        wHigh = actHigh/(actLow + actHigh);
        wLow  = actLow/(actLow + actHigh);

        % Total chemical ignition delay in [s]
        tauChem = (wLow*tauLow + wHigh*tauHigh)*1e-3;
    otherwise
        tauChem = 0;
end

end

function [y] = evalArrheniusEqn(theta,p,xiO2,ArrPar)

y =...
    ArrPar.k*...
    (p/ArrPar.pRef)^ArrPar.nP*...
    (xiO2/ArrPar.xiO2Ref)^ArrPar.nO2*...
    exp(ArrPar.thetaEa/theta);

end