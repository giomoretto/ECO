function [IntVar] = parOPDef(parModel,varargin)
% This script calculates the total Intake Mass via Signals from the
% HFM-sensor, information is expected in the InputVar-struct, else MDot
% gets calculated via a fit from OPRef

if nargin == 2
    meas = varargin{1};
    NumSym             = 'Sym';
    
    InputVar.engSpd    = meas.Ne; % 1/s; 2000 1/min
    InputVar.pIM       = meas.p_im; %Pa 150000 Pa
    InputVar.thetaIM   = meas.T_im; %K 55°C
    InputVar.xiBgIM    = meas.x_bg; %
    InputVar.pRail     = meas.p_rail; %Pa 1000bar
    % assumption on exhaust and residual gas
    InputVar.xiRes     = 0.05; %
    InputVar.lambdaExh = 2;
    % NOT USED
    InputVar.pExh      = InputVar.pIM + 0.3e5;%Pa
    InputVar.thetaExh  = 430+273.15; %°K
    % valve timings
    InputVar.ivc       = -172; % degCAaTDC
    InputVar.evo       = 155; % degCAaTDC
    InputVar.evc       = 345; % degCAaTDC
else
    
    NumSym             = 'Num';
    %% Define Engine Operation Parameters
    InputVar.engSpd    = 2000/60; % 1/s; 2000 1/min
    InputVar.pIM       = (1.2)*1e5; %Pa 1.5bar
    InputVar.thetaIM   = 318.15; %K 45°C
    InputVar.xiBgIM    = 0.1; %
    InputVar.pRail     = 1000e5; %Pa 1000bar
    % assumption on exhaust and residual gas
    InputVar.xiRes     = 0.05; %
    InputVar.lambdaExh = 2;
    % NOT USED
    InputVar.pExh      = InputVar.pIM + 0.3e5;%Pa
    InputVar.thetaExh  = 430+273.15; %°K
    % valve timings
    InputVar.ivc       = -172; % degCAaTDC
    InputVar.evo       = 155; % degCAaTDC
    InputVar.evc       = 345; % degCAaTDC
    
end

if nargin == 3
    meas = varargin{1};
    NumSym = varargin{2};
    
    InputVar.engSpd    = meas.Ne; % 1/s; 2000 1/min
    InputVar.pIM       = meas.p_im; %Pa 150000 Pa
    InputVar.thetaIM   = meas.T_im; %K 55°C
    InputVar.xiBgIM    = meas.x_bg; %
    InputVar.pRail     = meas.p_rail; %Pa 1000bar
    % assumption on exhaust and residual gas
    InputVar.xiRes     = 0.05; %
    InputVar.lambdaExh = 2;
    % NOT USED
    InputVar.pExh      = InputVar.pIM + 0.3e5;%Pa
    InputVar.thetaExh  = 430+273.15; %°K
    % valve timings
    InputVar.ivc       = -172; % degCAaTDC
    InputVar.evo       = 155; % degCAaTDC
    InputVar.evc       = 345; % degCAaTDC
end

%% Calculation with stationarity assumption
% Assumptions: complete combustion, i.e. only burnt gas and air are
% considered as components in exhaust gas
 % Volume at IVC
   [VInt,~]    = CylVol(InputVar.ivc,parModel.eng);
 % In-cylinder pressure at IVC
   pIVC = gasDynamicEffects(InputVar.pIM);
   
 % obtain fresh charge into cylinder and total mass -----------------------  
 % Fresh Charge into cylinder  
   mBeta = calcMBeta(InputVar,VInt);
 % account for resiudal gas fraction in total mass
   mCylTot  = mBeta / (1-InputVar.xiRes);

 % obtain residual mass in cylinder ---------------------------------------  
 % residual mass
   mResCyl = mCylTot*InputVar.xiRes;
 % Exhaust Gas Burnt Gas & Air Ratio Determination
   xiExhBg  = 1/InputVar.lambdaExh;
   xiExhAir = 1 - xiExhBg;
   
 % calculate composition of gas at Intake Valve Closing (IVC)--------------
 % Fresh Air Mass
   mBetaAir  = mBeta*(1-InputVar.xiBgIM);
 % Fresh Burnt Gas
   mBetaBg = mBeta-mBetaAir;
 % total fresh air
   mCylAir = mBetaAir + mResCyl*xiExhAir;
   mCylBg  = mBetaBg + mResCyl*xiExhBg;
 % In-cylinder species masses
   mCylN2   = parModel.thermo.xi.air.N2*mCylAir + parModel.thermo.xi.bg.N2*mCylBg;
   mCylO2   = parModel.thermo.xi.air.O2*mCylAir;
   mCylCxHy = 0;
   mCylCO2  = parModel.thermo.xi.bg.CO2*mCylBg;
   mCylH2O  = mCylTot - mCylN2 - mCylO2 - mCylCxHy - mCylCO2;

   xiCylAir  = mCylAir/mCylTot;
   xiCylBg   = 1 - xiCylAir;
   
   xiCylN2   = mCylN2/mCylTot;
   xiCylO2   = mCylO2/mCylTot;
   xiCylCxHy = mCylCxHy/mCylTot;
   xiCylCO2  = mCylCO2/mCylTot;
   xiCylH2O  = 1 - xiCylN2 - xiCylO2 - xiCylCxHy - xiCylCO2;

%% Initialization of output structure
   IntVar          = struct();
   IntVar.mCylTot  = mCylTot;
   IntVar.pInt     = pIVC;
   IntVar.VInt     = VInt;
   IntVar.xiAir    = xiCylAir;
   IntVar.xiBg     = xiCylBg;
   IntVar.xiO2     = xiCylO2;
   IntVar.xiN2     = xiCylN2;
   IntVar.xiCxHy   = xiCylCxHy;
   IntVar.xiCO2    = xiCylCO2;
   IntVar.xiH2O    = xiCylH2O;
   IntVar.ca_ivc   = -172; % degCAaTDC
   IntVar.ca_evo   = 155; % degCAaTDC
   IntVar.ca_evc   = 345; % degCAaTDC
   IntVar.engSpd   = InputVar.engSpd;
   IntVar.pRail    = InputVar.pRail;
   IntVar.pIM      = InputVar.pIM;
   IntVar.thetaIM  = InputVar.thetaIM;
   IntVar.xiBgIM   = InputVar.xiBgIM;
   [IntVar.kappaIvc,IntVar.specRIvc,IntVar.thetaIvc] = ...
       StaticCylinderConditions(0,pIVC,VInt,parModel,IntVar,NumSym);

end


function [mBeta,mDotBeta] = calcMBeta(InputVar,VInt)
pIM = InputVar.pIM;
thetaIM = InputVar.thetaIM;
Ne = InputVar.engSpd;

%% volumetric efficiency
% mean value for different pIVO/pIM
lambda_lp = 0.96;
% normalized for 45degC (318K) to one, measured for up to thetaIM = 90 degC
lambda_lTheta = (0.594 + 9.48e-4*thetaIM)/0.8955;
% measured only for rpm 2000
lambda_lOmega = 0.8955;
% combinded volumetric efficiencies
lambda_l = lambda_lp*lambda_lTheta*lambda_lOmega;

%% calculate mdot Beta
% ideal gas constant
R = 287;
% density of gas in intake manifold
rho = pIM/(R*thetaIM);
% frequency
freq = Ne/2;
% result
mBeta = lambda_l*rho*VInt;
% flow rate
mDotBeta = mBeta * freq;


end

function pIVC = gasDynamicEffects(pIM)

% everything in [Pa]
% pIVC = -3933 + 1.1003*pIM;
pIVC = -7122.9 + 1.1517*pIM;



end

