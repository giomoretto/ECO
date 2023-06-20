function [V,dVdeg,stroke] = CylVol(ca,eng)
% INPUTS:
%   - ca:  in [degCA] with TDC at ca = 0
%
% OUTPUTS:  
%   - V:   cylinder volume in [m^3]
%   - dV:  cylinder volume derivative in [m^3/deg]

% Define ca in radiant
  carad = ca*pi/180;
% Cylinder Volume
  V = eng.volDis /(eng.epsilon-1) + ...
      eng.volDis*1/2*...
      ( (eng.lengthConRod + eng.radCrank)/eng.radCrank - ...
        sqrt(eng.lengthConRod^2/eng.radCrank^2-sin(carad).^2) - ...
        cos(carad)...
      );
  
% Cylinder Volume Derivative
  dVrad = eng.volDis*1/2*...
                     ( sin(carad) + ...
                       (eng.lengthConRod^2/eng.radCrank^2-sin(carad).^2).^(-1/2) .* sin(carad).* cos(carad)...
                      );
  dVdeg = dVrad*pi/180;
  
% Piston Stroke Length, with rod ratio and % Stroke function f (Pischinger) 
  rodRatio = eng.radCrank/eng.lengthConRod;
  f = 1 + ...
      1/rodRatio - ...
      cos(carad) - ...
      1/rodRatio.*sqrt(1 - rodRatio^2*sin(carad).^2);
  stroke = eng.radCrank*f;

end