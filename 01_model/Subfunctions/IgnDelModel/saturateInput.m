function [y] = saturateInput(u,uMin,uMax,varargin)

saturationFun = 'smooth';

if nargin > 3
    saturationFun = varargin{1};
else
end

y = u;

switch saturationFun
    case 'exact'
        idxMin = (u <= uMin);
        idxMax = (u >= uMax);
        
        y(idxMin) = uMin;
        y(idxMax) = uMax;
        
    case 'smooth'
        epsilon = 1e-6;
        
        h = (uMax - uMin)/2;
        yMean = (uMin + uMax)/2;
        
        y(1:end) =...
            h/2*(sqrt(epsilon + ((u - yMean)/h + 1).^2) -...
                 sqrt(epsilon + ((u - yMean)/h - 1).^2)) + yMean;
        
    otherwise
end

end