function fxInit = createAcadosSimulation_xInit(fuxdot,parOpt,parModel)

%% Prepare optimization
% Parameters
  nInputs    = parModel.nInputs;
  nStates    = parModel.nStates;
  CAs        = parOpt.deltaPhi;         % represents sampling time
  CAsAcados  = parOpt.deltaPhiAcados;   % represents sampling time
  nInt       = parOpt.nInt(1);
  nIntAcados = parOpt.nInt(end);
  Int        = parOpt.Int;

  %% Create Integrator Functions for Optimization Range
% Initialization
  import casadi.* 

% Definition of possible integrators
  xStart    = MX.sym('xStart',nStates+nInputs+1,1); % Initial States
  pStart    = MX.sym('pStart',5,1);                 % Operating Point
  CAInt     = CAs/nInt;                             % Step size of each RK4 interval

% Define Explicit Integration allowing several steps (Euler Forward and
% Runge Kutta 4)
  xEndRK    = xStart;
  xEndEF    = xStart;
  for l = 1:nInt
      k1x    = fuxdot(xEndRK, pStart);
      k2x    = fuxdot(xEndRK + CAInt/2 * k1x, pStart);
      k3x    = fuxdot(xEndRK + CAInt/2 * k2x, pStart);
      k4x    = fuxdot(xEndRK + CAInt * k3x, pStart);
      xEndRK = xEndRK + CAInt/6 * (k1x + 2*k2x + 2*k3x + k4x);
      xEndEF = xEndEF + CAInt * k1x;
  end % for

% Create CasADi functions
  FIntRK    = Function('FIntRK',{xStart,pStart}, {xEndRK});
  FIntEF    = Function('FIntEF',{xStart,pStart}, {xEndEF});
  
  %% Create Integrator Functions for Simulation until EVO
% Definition of possible integrators
  xStart    = MX.sym('xStart',nStates+nInputs+1,1); % Initial States
  pStart    = MX.sym('pStart',5,1);                 % Operating Point
  CAInt     = CAsAcados/nIntAcados;                 % Step size of each RK4 interval

% Define Explicit Integration allowing several steps (Euler Forward and
% Runge Kutta 4)
  xEndRK    = xStart;
  xEndEF    = xStart;
  for l = 1:nIntAcados
      k1x    = fuxdot(xEndRK, pStart);
      k2x    = fuxdot(xEndRK + CAInt/2 * k1x, pStart);
      k3x    = fuxdot(xEndRK + CAInt/2 * k2x, pStart);
      k4x    = fuxdot(xEndRK + CAInt * k3x, pStart);
      xEndRK = xEndRK + CAInt/6 * (k1x + 2*k2x + 2*k3x + k4x);
      xEndEF = xEndEF + CAInt * k1x;
  end % for

% Create CasADi functions
  FIntRKAcados  = Function('FIntRKAcados',{xStart,pStart}, {xEndRK});
  FIntEFAcados  = Function('FIntEFAcados',{xStart,pStart}, {xEndEF});
  
%% Generate Integrating Function
% Create model states and inputs
  ux        = MX.sym('ux',nStates+nInputs+1,1);
  p         = MX.sym('p',5,1);

  % Run model integrator for optimization range
  caNum     = parOpt.ca(1:end-numel(parOpt.caAcados));
  xOut      = ux;
  for k = 1:numel(caNum)-1
      switch Int
          case 'EulerFW'
              xOut = [xOut, FIntEF(xOut(:,end),p)];
          case 'RK4'
              xOut = [xOut, FIntRK(xOut(:,end),p)];
          otherwise
              xOut = [xOut, FIntRK(xOut(:,end),p)];
      end
  end
  % Run model integrator for simulation until EVO
  caNum     = parOpt.ca(k+1:end);
  for k = 1:numel(caNum)-1
      switch Int
          case 'EulerFW'
              xOut = [xOut, FIntEFAcados(xOut(:,end),p)];
          case 'RK4'
              xOut = [xOut, FIntRKAcados(xOut(:,end),p)];
          otherwise
              xOut = [xOut, FIntRKAcados(xOut(:,end),p)];
      end
  end
  
% Reshape ux_trajectory
  xOut        = xOut(:);

% Define Model Functions
  fxInit      = Function('fxInit',{ux,p},{xOut},{'ux','p'},{'xOut'});
  fxInit      = fxInit.expand();

end
