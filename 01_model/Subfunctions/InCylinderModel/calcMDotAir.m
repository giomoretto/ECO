function mDotAir  = calcMDotAir(pCylIvc,xEGR)
% linear fit on Data from RefOP DataSet to calculate Massflow of Air into
% cylinder in dependence of pCyl [Pa] at IVC and xEGR [-]

%Data is fitted between 1.1e5<pCyl<2e5, 0<xEGR<0.5 

p1 = 0.01391;
p2 = 1.625e-07;
p3 = -0.03306;

% if xEGR < 0 
%     warning('EGR-ratio is smaller than 0')
% end

mDotAir = p1 + p2*pCylIvc + p3*xEGR;


end

