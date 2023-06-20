function [gamma] = evalCombWeightingFun(tau,tauMin,tauMax,varargin)

weightingFun = 'smooth';

if nargin > 3
    weightingFun = varargin{1};
else
end

switch weightingFun
    case 'non-smooth'
        if (tau < tauMin)
            gamma = 0;
        elseif (tau > tauMax)
            gamma = 1;
        else
            deltaTau = tauMax - tauMin;

            gamma = 1/deltaTau*tau - tauMin/deltaTau;
        end
        
    case 'smooth'
        slope = 1/(tauMax - tauMin);
        shift = (tauMax - tauMin)/2 + tauMin;

        gamma = 1/(1 + exp(-4*slope*(tau - shift)));
        
    otherwise
        slope = 1/(tauMax - tauMin);
        shift = (tauMax - tauMin)/2 + tauMin;

        gamma = 1/(1 + exp(-4*slope*(tau - shift)));
        
end

end

