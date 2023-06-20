function [tauIgn] = ignDelModel(pCyl,thetaCyl,xiO2,parModel)

% Inputs:   pCyl:       Actual Cylinder pressure [Pa]
%           thetaCyl:   Actual Cylinder Temperature [°K]
%           Par:        Parameter Struct
% Output:   tauIgn:     Ignition Delay [s]

    % Ignition delay time in [s]
    % Chemical ignition delay time in [s]
    switch parModel.Opts.IgnitionDelayModel
                    
        case 'Joerg' %*1.25 for xBG = 30%
            tauIgn_chem  = ignitionDelayJoerg(pCyl,thetaCyl,xiO2)*1e-3;
        case 1 % Simulink does not allow character arrays or strings
            tauIgn_chem  = ignitionDelayJoerg(pCyl,thetaCyl,xiO2)*1e-3;
        otherwise
            error('Wrong specification for ignition delay model.');
    end

    % Physical ignition delay time in [s]
    tauIgn_phys = parModel.ign.phys.const;

    % Constant bound
    tauIgnMin = parModel.ign.tauMin;
    tauIgnMax = parModel.ign.tauMax;

    % Total ignition delay time in [s]
    tauIgn =...
        saturateInput((tauIgn_chem + tauIgn_phys),...
                      tauIgnMin,tauIgnMax);
                
end

