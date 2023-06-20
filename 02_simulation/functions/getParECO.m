function parECO = getParECO(parOpt,parSim,parModel)
% this is a quite messy function which extracts all values needed for the
% simulink simulation and combines them into one single struct. In future,
% this struct should be a little cleaner and better organized

% Initialization
  parECO           = struct;
  
% Save other parameter structures needed for simulink visualizations
  parECO.parModel  = parModel;
  parECO.parSim    = parSim;
  
% Add numeric values needed for optimization
  nInj             = parModel.nInputs/2;
  SOENames         = cell(1,nInj); SOENames(1,:) = {'SOE'};
  DOENames         = cell(1,nInj); DOENames(1,:) = {'DOE'};
  uNames           = [SOENames,DOENames];

  parECO.u0        = parOpt.u0;
  parECO.Reference = parOpt.Reference;
  parECO.N         = numel(parOpt.ca) - 1;
  parECO.NAcados   = numel(parOpt.caAcados) - 1;
  parECO.T_ECO     = 0.06;
  parECO.nIter     = 1;
  parECO.scaleU    = [repmat(parOpt.Scale(end-1),parModel.nInputs/2,1); ...
                      repmat(parOpt.Scale(end),parModel.nInputs/2,1)];
  parECO.offsU     = [repmat(parOpt.Offs(end-1),parModel.nInputs/2,1); ...
                      repmat(parOpt.Offs(end),parModel.nInputs/2,1)];
  parECO.scaleX    = parOpt.Scale(1:parModel.nStates)';
  parECO.offsX     = parOpt.Offs(1:parModel.nStates)';
  parECO.scaleCA   = parOpt.Scale(end-1);
  parECO.offsCA    = parOpt.Offs(end-1);
  parECO.ca        = parOpt.ca;
  parECO.ca_ivc    = parSim.OP.ca_ivc;
  parECO.lbu       = scaleUnscale(parOpt.uMin,uNames,parOpt,'scale',false);
  parECO.ubu       = scaleUnscale(parOpt.uMax,uNames,parOpt,'scale',false);
  parECO.idxCoC    = find(parOpt.ca >= parOpt.Reference.CoCmax,1,'first');
  
% Add numeric values needed for simulation
% Calculate x0 and initial trajectory with simulation
  parECO.ux0       = calcInit(parOpt,parSim,parModel);
  parECO.ca_x0Opt  = parSim.OP.ca_ivc : 4 : parOpt.OptimizationRange(1);
  parECO.dca_x0Opt = parECO.ca_x0Opt(2)-parECO.ca_x0Opt(1);
  
% Simulate plant to generate measurements
  parECO.ca_sim    = parOpt.OptimizationRange(1) : 0.5 : parSim.OP.ca_evo;
  parECO.dca_sim   = parECO.ca_sim(2)-parECO.ca_sim(1);
  
  parECO.dca_opt       = parOpt.deltaPhi;
  parECO.dca_optAcados = parOpt.deltaPhiAcados;

% Add numeric values needed for model
  parECO.nInputs   = parModel.nInputs;
  parECO.nStates   = parModel.nStates;
  parECO.nOutputs  = parModel.nOutputs;
  parECO.nConstr   = 3;
  
end

function ux0Opt = calcInit(parOpt,parSim,parModel)
  % Prepare Simulation
    caPreOpt    = (parSim.OP.ca_ivc:parSim.Opts.deltaPhi:...
                        parOpt.OptimizationRange(1));
    u0PreSim    = repmat(parOpt.u0,1,numel(caPreOpt));
    x0PreOpt    = [parSim.OP.pInt,0,0]';            % [p QComb IMEP]
    if parOpt.enNOx
        x0PreOpt = [x0PreOpt;parSim.OP.thetaIvc;0]; % [p QComb IMEP Tuz cNOx]
    end

  % Run simulation
    simPreOpt   = CompleteSimulation(caPreOpt,x0PreOpt,u0PreSim,parSim,parModel);

  % Get x0 at Start of Opt Range and build ux0
    x0Opt       = simPreOpt.x(:,end);
    ux0Opt      = [parOpt.u0; x0Opt; parOpt.OptimizationRange(1)];
    
end