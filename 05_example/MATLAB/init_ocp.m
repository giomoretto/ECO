clear mex; %#ok<CLMEX>
buildFolder       = what('06_build');
acadosBuildFolder = fullfile(buildFolder.path,'acadosBuild');
if isfolder(acadosBuildFolder)
    rmpath(genpath(acadosBuildFolder));
    rmdir(acadosBuildFolder,'s');
    pause(0.5);
end

%% Parameter Definition
% Define Optimization Parameters
parOpt.Reference.IMEP    = 6e5;
parOpt.Reference.pMax    = 150e5;
parOpt.Reference.dpMax   = 4e5;
parOpt.Reference.Tmin    = 0+273;
parOpt.Reference.cNOx    = 1e4;
parOpt.Reference.CoCmax  = 20;
parOpt.Reference.PhiMax  = 1/1.3;
parOpt.deltaPhi          = 0.5;
parOpt.deltaPhiAcados    = 3;           % for acados framework, precision from 
                                        % end of optimization range until EVO
parOpt.OptimizationRange = [-16 26];    % °CA ATDC
parOpt.Int               = 'RK4';       % RK4 / EulerFW / EulerBW / EulerMix
parOpt.nInt              = 1;
parOpt.nIntAcados        = 1;
parOpt.SQP.StepSize      = 1;
parOpt.SQP.nSQPmax       = 100;
parOpt.enNOx             = true;
parOpt.ca                = (parOpt.OptimizationRange(1):...
    parOpt.deltaPhi:...
    parOpt.OptimizationRange(2));
parOpt.u0                = [-10.3;-5.3;...
    167;372];
for i = 1:numel(parOpt.u0)/2
    parOpt.uMin(i) = parOpt.OptimizationRange(1);
    parOpt.uMax(i) = 10; 
    parOpt.uMin(i+numel(parOpt.u0)/2) = 80;
    parOpt.uMax(i+numel(parOpt.u0)/2) = 800;
end
parOpt.uMin = parOpt.uMin(:);
parOpt.uMax = parOpt.uMax(:);

parOpt.xMin              = [0;0;-2e6];
parOpt.xMax              = [parOpt.Reference.pMax;2e3;2e6];
parOpt.dt_Inj            = 400; % min. distance between injections

% Variable Scaling
% Def: scaled = (unscaled-Offs)/Scale -> unscaled = scaled*Scale+Offs
parOpt.ScaleOffsVars = {'pCyl','QComb','IMEP','Theta','NO','NOppm','SOE','DOE'};
parOpt.Offs          = [ 0,     0,      0,     0,      0,    0,     -20,  100];
parOpt.Scale         = [ 100e5, 1000,  10e5,   1e3,    1e-6, 1e3,    70,  700];

% load parameters
parModel = ModelParameters();
parModel.nInputs = numel(parOpt.u0);

% Include NOx model
if parOpt.enNOx
    parOpt.xMin       = [parOpt.xMin;0;0];
    parOpt.xMax       = [parOpt.xMax;1e4;0.1];
    parModel.nStates  = 5;
    parModel.nOutputs = 12;
    parModel.NOxmodel = 'xBZ'; % xBG xBZ
end

% Define Engine Simulation Parameters
parSim.OP            = parOPDef(parModel);
parSim.Opts.deltaPhi = 0.5;
parSim.Opts.Int      = 'RK4'; % RK4 / EulerFW
parSim.Opts.nInt     = 1;
parSim.Opts.ca       = (parSim.OP.ca_ivc:parSim.Opts.deltaPhi:...
    parSim.OP.ca_evo);

%% Create Casadi Functions and Run Solver
% Extend optimization range until EVO for acados framework
parOpt.caAcados = parOpt.ca(end)+parOpt.deltaPhiAcados : ...
    parOpt.deltaPhiAcados : parSim.OP.ca_evo;
parOpt.nInt     = [repmat(parOpt.nInt,1,numel(parOpt.ca)-1), ...
    repmat(parOpt.nIntAcados,size(parOpt.caAcados))];
parOpt.ca       = [parOpt.ca, parOpt.ca(end)+parOpt.deltaPhiAcados : ...
    parOpt.deltaPhiAcados : parSim.OP.ca_evo];

% Create ocp and solve
ocp = createAcadosFunctions_InjOpt(parOpt,parModel);
runSQPacados_InjOpt(ocp,parModel,parSim,parOpt);
