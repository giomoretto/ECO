function fsim = createAcadosSimulation(fuxdot,fy,parOpt,parModel)

%% Prepare optimization
% Parameters
  nInputs   = parModel.nInputs;
  nStates   = parModel.nStates;
  nInj      = nInputs/2;
  nInt      = parOpt.nInt(1);
  Int       = parOpt.Int;
  if parOpt.enNOx
      xNames  = {'pCyl','QComb','IMEP','Theta','NO'};
      yNames  = {'Theta','kappa','dQWall','specR','mFuelPrep','tauIgn','xBG','Theta','NOppm'};
  else
      xNames  = {'pCyl','QComb','IMEP'};
      yNames  = {'Theta','kappa','dQWall','specR','mFuelPrep','tauIgn'};
  end
  SOENames  = cell(1,nInj); SOENames(1,:) = {'SOE'};
  DOENames  = cell(1,nInj); DOENames(1,:) = {'DOE'};
  uxNames   = [SOENames, DOENames, xNames, {'SOE'}];

%% create Integrator Functions
% Initialization
  import casadi.* 

% Definition of possible integrators
  xStart    = MX.sym('xStart',nStates+nInputs+1,1); % Initial States
  pStart    = MX.sym('pStart',6,1);                 % Operating Point
  dCA       = MX.sym('dCA',1,1);                    % Crank-angle Resolution
  CAInt     = dCA/nInt;                             % Step size of each RK4 interval

% Define Explicit Integration allowing several steps (Euler Forward and
% Runge Kutta 4)
  xEndRK    = scaleUnscale(xStart,uxNames,parOpt,'scale',false);
  xEndEF    = scaleUnscale(xStart,uxNames,parOpt,'scale',false);
  for l = 1:nInt
      k1x    = fuxdot(xEndRK, pStart);
      k2x    = fuxdot(xEndRK + CAInt/2 * k1x, pStart);
      k3x    = fuxdot(xEndRK + CAInt/2 * k2x, pStart);
      k4x    = fuxdot(xEndRK + CAInt * k3x, pStart);
      xEndRK = xEndRK + CAInt/6 * (k1x + 2*k2x + 2*k3x + k4x);
      xEndEF = xEndEF + CAInt * k1x;
  end

% Create CasADi functions
  xEndRK    = scaleUnscale(xEndRK,uxNames,parOpt,'unscale',false);
  xEndEF    = scaleUnscale(xEndEF,uxNames,parOpt,'unscale',false);
  
  FIntRK    = Function('FIntRK',{xStart,pStart,dCA}, {xEndRK});
  FIntEF    = Function('FIntEF',{xStart,pStart,dCA}, {xEndEF});
  
%% Generate Integrating Function
% Create model states and inputs
  ux        = MX.sym('ux',nStates+nInputs+1,1);
  p         = MX.sym('p',6,1);
  dca       = MX.sym('dca',1,1);

% Run model integrator
  switch Int
      case 'EulerFW'
          xOut = FIntEF(ux,p,dca);
      case 'RK4'
          xOut = FIntRK(ux,p,dca);
      otherwise
          xOut = FIntRK(ux,p,dca);
  end
  
% Create and unscale state derivative and model outputs
  xOutScale = scaleUnscale(xOut,uxNames,parOpt,'scale',false);
  xdotScale = fuxdot(xOutScale,p);
  yOutScale = fy(xOutScale,p);
  
  xdot      = scaleUnscale(xdotScale,uxNames,parOpt,'unscale',true);
  yOut      = scaleUnscale(yOutScale,yNames,parOpt,'unscale',false);
  
% Define Model Functions
  fsim      = Function('fsim',{ux,p,dca},{xOut,xdot,yOut},...
                {'ux','p','dca'},{'xOut','xdot','yOut'});
  fsim      = fsim.expand();

end
