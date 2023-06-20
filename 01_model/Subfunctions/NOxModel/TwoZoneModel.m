function [dThetaUZ,thetaBZ] = TwoZoneModel(QComb,thetaUZ0,dpCyl,pCyl,dVCyl,VCyl,...
        TCyl,kappa,mTot,parOP,parModel,xBG,xBZ)

% global test
%TWOZONEMODEL calculates the burned gas temperature which is needed for the
%             NOx-model
%   This model assumes two zones, namely the unburned and the burned zone.
%   It calculates the unburned zone temperature under the assumption of
%   isentropic compression/expansion and constant gas properties.
%   Furthermore, a logarithmic barrier function is used to smoothly
%   approximate the unburned zone temperature throughout the whole high
%   pressure cycle. The burned zone temperature is then derived via the
%   average cylinder temperature.

%% Derive dThetaUZ and thetaBZ
 % Limit ThetaUZ to TCyl at current shooting node
   thetaUZ0 = saturateInput(thetaUZ0,0,TCyl);

%    uMin     = 5;
%    uMax     = 10;
%    uComb    = (saturateInput(QComb,uMin,uMax) - uMin)/(uMax-uMin);
   dTheta   = TCyl*(dpCyl/pCyl + dVCyl/VCyl);
   zetaComb = QComb/(mTot*parModel.thermo.fuel.lowHeatVal);


   Qmin = 0.01;
   Qmax = 0.1;
   uComb     = (saturateInput(zetaComb,Qmin,Qmax)-Qmin)/(Qmax-Qmin);

 % unburned zone temperature
   dThetaUZ = 1/pCyl*dpCyl*thetaUZ0*(kappa-1)/kappa;
   dThetaUZ = (1-uComb)*dTheta + uComb*dThetaUZ;
   
 % limit ThetaUZ to TCyl
%    uTheta       = saturateInput(100*(TCyl+0.25*dTheta-(thetaUZ0+0.25*dThetaUZ)),0,1);
%    dThetaUZ     = (1-uTheta)*dTheta + uTheta*dThetaUZ;
%    thetaUZ0     = (1-u)*TCyl + u*thetaUZ0;

 % Limit future ThetaUZ to TCyl using euler forward integration
%    deltaPhi = 0.05;
%    dThetaUZ = saturateInput(dThetaUZ,-100,(TCyl+deltaPhi*dTheta-thetaUZ0)/deltaPhi);
   
 % burned zone temperature
%    mAir         = mTot*parModel.thermo.fuel.airFuelRatioSt;
%    mg_stoich    = mAir/parOP.xiAir;
%    zeta_comb    = QComb/(mTot*parModel.thermo.fuel.lowHeatVal);
%    xBZ          = (1-uComb)*1e3 + uComb*zeta_comb * mg_stoich/parOP.mCylTot;
%    thetaBZ      = (TCyl-thetaUZ0)/xBZ + thetaUZ0;
%    thetaBZ  = (TCyl-thetaUZ0)/xBG + thetaUZ0;
   xBZsat = saturateInput(xBZ,0.01,1);
   thetaBZ  = uComb*(TCyl-thetaUZ0)/xBZsat + thetaUZ0;

%    test.xBZ(test.i) = xBZ;
%    test.xBG(test.i) = xBG;
%    test.xBZsat(test.i) = xBZsat;
end
