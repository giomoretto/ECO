function [dNOdphi,NOppm,NOeq] = NOxModel(QComb,thetaUZ,thetaBZ,NO,mTot,parOP,parModel,xBG,VCyl,dVCyl,pCyl,dpCyl,thetaCyl)
%NOXMODEL calculates the NOx concentration
%   Output is mole franction in the burned gas. In order to be compared with
%   measured NOx concentrations the output has to be multiplied by 1e6*xBG

% global test

switch parModel.NOxmodel
    case 'xBG'
%% Derive dNOdphi
   d       = parModel.NoParams;
   lambda  = parOP.mCylTot*parOP.xiAir/mTot/parModel.thermo.fuel.airFuelRatioSt;
   lambda  = saturateInput(lambda,1,10);
   
 % Time passed per °CA
   dTdPhi  = 1/(parOP.engSpd*360);

% NO-equilibrium concentration
   NOeq    = NOequilibrium(thetaBZ, 1); % Mole fraction in BG
%    NOeq = xBG*NOeq; % Mole fraction in whole gas

 % Arrhenius
   R1      = d(1).*exp(d(2).*lambda).*exp(d(3)./thetaBZ);
   Q       = d(4).*exp(d(5)./thetaBZ);
%    R1      = d(1).*exp(d(2)./thetaBZ);
%    Q       = d(3).*exp(d(4)./thetaBZ);
   
 % Formula from David Ueberschär, who adapted the formula from Eriksson (2016)
   uMin     = 1795;
   uMax     = 1800;
   uNOx     = (saturateInput(thetaBZ,uMin,uMax) - uMin)/(uMax-uMin);

%    u       = (saturateInput(thetaBZ,1795,1800)-1795)/5;
   dNOdt   = uNOx * R1.*(1-(NO/NOeq).^2)/(1+NO/NOeq * Q);
   dNOdphi = dNOdt*dTdPhi; % Mole fraction in BG
%    dNOdphi = xBG*dNOdphi; % Mole fraction in whole gas

% dNOdphi is assumed to only take postivie values, to prevent the solver to
% reach negative NO values.
% eps        = 1e-15;
% dNOdphi    = 1/2*dNOdphi + 1/2*sqrt(eps + dNOdphi^2);
eps        = 1e-8;
dNOdphi    = 1/2*dNOdphi + 1/2*sqrt(eps^2 + dNOdphi^2);%-0.5*eps;

%  % Include volume change
% %    mIVC    = parOP.mCylTot;
% %    mUZ     = mIVC - mIVC*mFuel/mTot/parModel.comb.etaComb;
% %    VUZ     = mUZ*287*thetaUZ/pCyl;
% %    VBZ     = VCyl-VUZ;
% %    dVUZ    = -dp*VUZ/kappa/pCyl;
% %    dVBZ    = dVCyl - dVUZ;
% %    dNOdphi = dNOdphi - NO/VBZ*dVBZ;
%    dNOdphi = dNOdphi - NO/VCyl*dVCyl;

 % Convert NO concentration to ppm in whole gas
   NOppm   = 1e6*(xBG-parOP.xiBg)*NO;
%    mAir         = mTot*parModel.thermo.fuel.airFuelRatioSt;
%    mg_stoich    = mAir/parOP.xiAir;
%    zeta_comb    = QComb/(mTot*parModel.thermo.fuel.lowHeatVal);
%    xBZ          = zeta_comb * mg_stoich/parOP.mCylTot;
%    NOppm        = 1e6*xBZ*NO;
%    NOppm        = 1e6*xBG*NO;
   
    case 'xBZ'
%% Derive dNOdphi
 % Load parameters
%    d       = parModel.NoParams;
      
 % Time passed per °CA
   dTdPhi  = 1/(parOP.engSpd*360);
   
 % Mole fraction to concentration in BG (mol/cm^3)
   mf2c    = 1e-6*pCyl/(parModel.thermo.gas.univR*thetaBZ);

 % Arrhenius
%    R1      = 4.6838e17*exp(-115423/thetaBZ);
%    R2      = 5.6706e16*exp(-80426/thetaBZ);
%    R3      = 2.4438e17*exp(-80310/thetaBZ);
%    Q       = 1.5925*exp(-35162/thetaBZ);
%    Q = R1/(R2+R3);

thetaBZsat = saturateInput(thetaBZ,1600,3300);
% thetaBZsat = thetaBZ;

% reaction rates
[k1f,k2f,k3f] = reacRates(thetaBZsat);
 % equilibrium concentrations   
X_N  = X_N_eq(thetaBZsat);
X_O  = X_O_eq(thetaBZsat);
X_N2 = X_N2_eq(thetaBZsat);
X_O2 = X_O2_eq(thetaBZsat);
X_OH = X_OH_eq(thetaBZsat);
X_NO = X_NO_eq(thetaBZsat);
NOeq = X_NO*mf2c;


R1 = k1f.*X_O.*X_N2.*mf2c.^2;
R2 = k2f.*X_N.*X_O2.*mf2c.^2;
R3 = k3f.*X_N.*X_OH.*mf2c.^2;

Q = R1./(R2+R3);
   
%    R1      = d(1).*exp(-d(3)./thetaBZ); % mol/(cm^3*s)
%    Q       = d(4).*exp(-d(5)./thetaBZ); % [-]
   
 % Turn on NOx model if Tbz > 1800K
%    uMin     = 1795;
%    uMax     = 1800;
%    uNOx     = (saturateInput(thetaBZ,uMin,uMax) - uMin)/(uMax-uMin);

 % Turn on NOx model if combustion has started
   zetaComb = QComb/(mTot*parModel.thermo.fuel.lowHeatVal);

   Qmin = 0.01;
   Qmax = 0.05;
   uNOx     = (saturateInput(zetaComb,Qmin,Qmax)-Qmin)/(Qmax-Qmin);
 % Formula from David Ueberschär, who adapted the formula from Eriksson (2016)
%    NOscale = saturateInput(NO/NOeq,0,1);
%    NOscale = saturateInput(-uNOx,1,NO/NOeq);
NOscale = NO/NOeq;
   dNOdt    = uNOx * 2*R1.*(1-(NOscale).^2)/(1+NOscale * Q);
%     dNOdt  = uNOx*2*R1;

 % Include volume change
   LHV      = parModel.thermo.fuel.lowHeatVal;
   sigma    = parModel.thermo.fuel.airFuelRatioSt;
   mUZ      = parOP.mCylTot - QComb/LHV*sigma*(1/(1-parOP.xiBg));
   VUZ      = mUZ*parOP.specRIvc*thetaUZ/pCyl;
   VBZ      = VCyl-VUZ;
   dVUZ     = -dpCyl*VUZ/parOP.kappaIvc/pCyl;
   dVBZ     = uNOx*(dVCyl - dVUZ);
   dNOdV    = uNOx*NO*dVBZ/VBZ;
   dNOdphi  = dNOdt*dTdPhi - 1*dNOdV;


 % dNOdphi is assumed to only take postivie values, to prevent the solver to
 % reach negative NO values.
%    eps      = 1e-8;
%    dNOdphi  = 1/2*dNOdphi + 1/2*sqrt(eps^2 + dNOdphi^2) - 0.5*eps;

 % Convert NO concentration to ppm in whole gas
%    maStoich = mTot*parModel.thermo.fuel.airFuelRatioSt;
%    mgStoich = maStoich/parOP.xiAir;
%    zetaComb = QComb/(mTot*parModel.thermo.fuel.lowHeatVal);
%    xBZ      = zetaComb * mgStoich/parOP.mCylTot;
%    NOppm    = 1e6*xBZ*NO;
%    NOppm    = 1e6*(xBG-parOP.xiBg)*NO;
   
   NOppm = 1e6*NO * (1e6*VBZ) / (pCyl*VCyl/(parModel.thermo.gas.univR*thetaCyl));
%    test.dVBZ(test.i) = dVBZ;
%    test.dVUZ(test.i) = dVUZ;
%    test.dNOchem(test.i) = dNOdt*dTdPhi;
%    test.dNOdV(test.i) = uNOx*NO*dVBZ/VBZ;
%    test.dNOdphi(test.i) = dNOdphi;
%    test.R1conc(test.i) = R1;
%    test.R2conc(test.i) = R2;
%    test.R3conc(test.i) = R3;
%    test.Qconc(test.i) = Q;
%    test.k1f(test.i) = k1f;
%    test.X_N2(test.i) = X_N2;
%    test.X_O(test.i) = X_O;
%    test.NOeq(test.i) = NOeq;
%    test.X_NO(test.i) = X_NO;
%    test.NOscale(test.i) = NOscale;
%    test.thetaBZ(test.i) = thetaBZ;
%    test.QComb(test.i) = QComb;
%    test.uNOx(test.i) = uNOx;
%    test.Vbz(test.i) = VBZ;
%    test.dVbz(test.i) = dVBZ;
%    test.zetaComb(test.i) = zetaComb;
%    test.mf2c(test.i) = mf2c;
%    test.nTot(test.i) = pCyl*VCyl/(parModel.thermo.gas.univR*thetaCyl);
%    test.i = test.i + 1;
end


    function y = NOequilibrium(T_burned, lambda)
        % calculates the Mole fraction of NOequilibrium at T_burned, lambda
        %        y0 = 0.0705;
        y0 = 0.14343;

        %        b1 = -156.2;
        %        b2 =  0.9122;
        %        b3 =  157.3;
        %        b4 =  0.9104;

        c1 =  21.57;
        c2 = -0.5528;
        c3 = -1.846;

        %        y  = y0 * (b1*exp(b2*lambda)+b3*exp(b4*lambda)) * exp(-c1*(T_burned/1000-c2)^c3);
        y  = y0 * exp(-c1.*(T_burned./1000-c2).^c3);

        %        % Cubic fit
        %        a1 = -1.6660e-12;
        %        a2 =  1.7847e-08;
        %        a3 = -4.2743e-05;
        %        a4 =  3.0544e-02;
        %
        %        y  = a1*T_burned.^3 + a2*T_burned.^2 + a3*T_burned + a4;
        %
        %        % Quadratic fit
        %        a1 =  5.4561e-09;
        %        a2 = -1.0810e-05;
        %        a3 = -6.2422e-03;
        %
        %        y  = a1*T_burned.^2 + a2*T_burned + a3;
        %
        %        % Linear fit
        %        % a1 =  2.4696e-05;
        %        % a2 = -6.3406e-02;
        %
        %        % y  = a1*T_burned + a2;


    end

    function y = NOequilibriumNew(Tbz,xBG)
        % Extract the mole fraction of NOequilibrium at Tbz and xBG from map

        % Limit input according to map
        xBG = saturateInput(xBG,0,1);

        % Interpolate NOeq from map
        %        NOeq_map = load('NOeq_map.mat');
        %        y = interp2(NOeq_map.Tbz,NOeq_map.xBG,NOeq_map.NOeq,Tbz,xBG,'cubic');

        % Estimate NOeq from quadratic fit
        a1 =  4.26602564102597e-12;
        a2 = -2.07548451548477e-08;
        a3 =  3.36296074758639e-05;
        a4 = -0.0180460099900153;
        y = a1*Tbz^3 + a2*Tbz^2 + a3*Tbz + a4;

    end

    function R1out = R1table(T_burned)

        temp = [1700, 1800:100:3000];
        R1tab = [0,   0.000000002179464,   0.000000020268221,   0.000000150078939,   0.000000913787326,   0.000004696472118,   0.000020818216974,   0.000081020245524,   0.000280967204499,   0.000879089808274,   0.002507764092364,   0.006581323249315,   0.016012694873531,   0.036363364439072];
        R1out = interp1(temp,R1tab,T_burned,'linear','extrap');

    end

    function R2out = R2table(T_burned)

        temp = [1700, 1800:100:3000];
        R2tab = [0,   0.000000000173726,   0.000000001840163,  0.000000015339860,   0.000000104054332,   0.000000590140700,  0.000002861223281,   0.000012077975187,   0.000045066530629,   0.000150530581582,   0.000454900260454,   0.001255001188166,   0.003185388470233,   0.007488705772492];
        R2out = interp1(temp,R2tab,T_burned,'linear','extrap');

    end

    function R3out = R3table(T_burned)

        temp = [1700, 1800:100:3000];
        R3tab = [0,   0.000000001100244,   0.000000011075358,   0.000000088029820,   0.000000571276296,   0.000003110726581,   0.000014534982550,   0.000059369039441,   0.000215283515130,   0.000702094134301,   0.002081928146211,   0.005665714362881,   0.014262865009634,   0.033441790260539];
        R3out = interp1(temp,R3tab,T_burned,'linear','extrap');

    end

    function [k1f,k2f,k3f] = reacRates(Tbz)

        k1f = 1.47e13.*Tbz.^0.3.*exp(-37885.8746./Tbz);
        k2f = 6.40e09.*Tbz.^1  .*exp(-3163.169283./Tbz);
        k3f = 3.80e13;

    end

    function Xout = X_H_eq(Tbz)

        coef = [-57.5467573413888,...
                  0.0403118313783464,...
                 -1.12950250476401e-05,...
                  1.17790762646850e-09];
        regr = [1, Tbz, Tbz.^2, Tbz.^3];

        Xout = exp(coef*regr');

    end

    function Xout = X_N_eq(Tbz)

        coef = [-93.5994804475774,...
            0.0635582925091314,...
            -1.78636544365512e-05,...
            1.86389452885879e-09];
        regr = [1, Tbz, Tbz.^2, Tbz.^3];

        Xout = exp(coef*regr');

    end

    function Xout = X_N2_eq(Tbz)

        coef = [-0.216167752677653,...
            -0.000141506777271925,...
            8.06135838153914e-08,...
            -1.54929244598273e-11];
        regr = [1, Tbz, Tbz.^2, Tbz.^3];

        Xout = exp(coef*regr');

    end

    function Xout = X_NO_eq(Tbz)

        coef = [-34.406106632447,...
                 0.0234949344561571,...
                -6.41287032472529e-06,...
                 6.35111932336585e-10];
        regr = [1, Tbz, Tbz.^2, Tbz.^3];

        Xout = exp(coef*regr');

    end

    function Xout = X_O_eq(Tbz)

        coef = [-63.6851353343809,...
                  0.0451581256769548,...
                 -1.25158990443408e-05,...
                  1.27610634104909e-09];
        regr = [1, Tbz, Tbz.^2, Tbz.^3];

        Xout = exp(coef*regr');

    end

    function Xout = X_O2_eq(Tbz)

        coef = [-33.6431381400466,...
                  0.0227657004453811,...
                 -6.04189495663701e-06,...
                  5.68577263026475e-10];
        regr = [1, Tbz, Tbz.^2, Tbz.^3];

        Xout = exp(coef*regr');
    end

    function Xout = X_OH_eq(Tbz)

        coef = [-38.9129885343657,...
                  0.0269739906905952,...
                 -7.40667321236856e-06,...
                  7.42212500799537e-10];
        regr = [1, Tbz, Tbz.^2, Tbz.^3];

        Xout = exp(coef*regr');
    end


end

