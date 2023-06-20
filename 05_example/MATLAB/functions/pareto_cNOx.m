function out = pareto_cNOx(ocp,parOpt,parModel,parSim)
% This function solves the ocp for a sequence of reference values for the
% maximum NOx Concentration (cNOx)

% Define reference values
  cNOx     = [1e4 1550 1200 900];

% Initialize model simulation
  caSim    = parSim.Opts.ca;
  x0Sim    = [parSim.OP.pInt,0,0]';
  % Include NOx model
  if parOpt.enNOx
      x0Sim = [x0Sim;parSim.OP.thetaIvc;0];
  end
  
% Initialize function outputs
  uOpt     = zeros(parModel.nInputs,numel(cNOx));
  y_eta    = zeros(1,numel(cNOx));
  y_cNOx   = zeros(1,numel(cNOx));
  idx_pCyl = find(caSim >= -10 & caSim <= 29);
  y_pCyl   = zeros(numel(idx_pCyl),numel(cNOx));
  status   = zeros(1,numel(cNOx));

% Loop over all reference values and solve ocp
  for k = 1:numel(cNOx)
      % Update reference value
      parOpt.Reference.cNOx = cNOx(k);
      % Solve ocp
      [uOpt(:,k),status(1,k)] = runSQPacados_InjOpt(ocp,parModel,parSim,parOpt);
      % Run simulation
      uSim        = repmat(uOpt(:,k),1,numel(caSim));
      simOpt      = CompleteSimulation(caSim,x0Sim,uSim,parSim,parModel);
      % Extract solution
      imep        = simOpt.x(3,end);
      mFuel       = simOpt.y(5,end);
      y_eta(k)    = 100*imep*parModel.eng.volDis/(mFuel*42.6*1e6);
      y_cNOx(k)   = simOpt.y(9,end);
      y_pCyl(:,k) = simOpt.x(1,idx_pCyl)' * 1e-5;
  end

% Create output struct
  out.name    = '$X_{\rm NO,max}$';
  out.uOpt    = uOpt;
  out.y       = y_cNOx;
  out.unit    = '[ppm]';
  out.eta     = y_eta;
  out.pCyl    = y_pCyl;
  out.ca_pCyl = simOpt.ca(idx_pCyl);
  out.status  = status;

end