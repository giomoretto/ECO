function x = checkValidity(x,ca,parModel,parOP)

% Check if NOx model is enabled
if numel(x) > 3
    enNOx = true;
else
    enNOx = false;
end

% Read states;
pCyl    = x(1);
QComb   = x(2);
imep    = x(3);
if enNOx
    Tuz = x(4);
    NOx = x(5);
end

% Limit states to a physical lower bound
% pCyl
eps     = 1e-3;
VCyl    = CylVol(ca,parModel.eng);
lb      = 0.8 * parOP.pInt*(parOP.VInt./VCyl).^1.4;
pCyl    = 1/2*(pCyl-lb) + 1/2*sqrt(eps^2 + (pCyl-lb)^2) + lb;

% max Qcomb value
xBZmax = 0.95;
QCombMax = (-parModel.thermo.fuel.lowHeatVal*parOP.mCylTot*xBZmax)/...
    (xBZmax - parModel.thermo.fuel.airFuelRatioSt/(1-parOP.xiBg)-1);
QComb = saturateInput(QComb,0,QCombMax);


if enNOx
    eps = 1e-10;
    lb  = 0;
    NOx = 1/2*(NOx-lb) + 1/2*sqrt(eps^2 + (NOx-lb)^2) + lb;% - eps/2;
end

% Create output vector
x       = [pCyl; QComb; imep];
if enNOx
    x   = [x; Tuz; NOx];
end

end