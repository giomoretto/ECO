%% Optimization Main Script
clear; clear mex; %#ok<CLMEX>
buildFolder       = what('06_build');
acadosBuildFolder = fullfile(buildFolder.path,'acadosBuild');
if isfolder(acadosBuildFolder)
    rmpath(genpath(acadosBuildFolder));
    rmdir(acadosBuildFolder,'s');
    pause(0.5);
end

warning('off','MATLAB:linkaxes:RequireDataAxes')
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

%% Parameter Definition
% Define Optimization Parameters
parOpt.Reference.IMEP    = 6e5;
parOpt.Reference.pMax    = 150e5;
parOpt.Reference.dpMax   = 4e5;
parOpt.Reference.Tmin    = 0+273;
parOpt.Reference.cNOx    = 1e4;
parOpt.Reference.CoCmax  = 20;
parOpt.Reference.PhiMax  = 1/1.3;
parOpt.Cost.WSlackdpMax  = 1e2;
parOpt.deltaPhi          = 0.5;
parOpt.deltaPhiAcados    = 3;			% for acados framework, precision from
										% end of optimization range until EVO
parOpt.OptimizationRange = [-16 26];	% °CA ATDC
parOpt.Solver            = 'Build';		% SQP / Build
parOpt.Int               = 'RK4';		% RK4 only integration method so far
parOpt.nInt              = 1;
parOpt.nIntAcados        = 1;
parOpt.SQP.StepSize      = 1;
parOpt.enSFunGen         = true;
% enable NOx model
parOpt.enNOx             = true;
parOpt.ca                = (parOpt.OptimizationRange(1):...
    parOpt.deltaPhi:...
    parOpt.OptimizationRange(2));
% initialize solution with two injections
parOpt.u0                = [-14.7; -9.6; ...
    80; 404];
% initialize solution with three injections, only possible within matlab
% parOpt.u0                = [-8; -5; -2; ...
%     150; 150; 150];
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
% check function scaleUnscale.m
parOpt.ScaleOffsVars = {'pCyl','QComb','IMEP','Theta','NO','NOppm','SOE','DOE'};
parOpt.Offs          = [ 0,     0,      0,     0,      0,    0,     -20,  100];
parOpt.Scale         = [ 100e5, 1000,  10e5,   1e3,    1e-6, 1e3,    70,  700];
% parOpt.Scale         = [ 100e5, 1000,  10e5,   1e3,   1e-7, 1e1,    7,  700];

% load Engine model and combustion model parameters
parModel = ModelParameters();  
parModel.nInputs = numel(parOpt.u0);

% Include NOx model
if parOpt.enNOx
    parOpt.xMin       = [parOpt.xMin;0;0];
    parOpt.xMax       = [parOpt.xMax;1e4;0.1];
    parModel.nStates  = 5;
    parModel.nOutputs = 12;
    parModel.NOxmodel = 'xBZ'; % xBG xBZ
else
    parModel.nStates = 3;
    parModel.nOutputs = 8;
end

% Define Engine Simulation Parameters
parSim.OP            = parOPDef(parModel);
parSim.Opts.deltaPhi = 0.5;
parSim.Opts.Int      = 'RK4'; % RK4 / EulerFW
parSim.Opts.nInt     = 1;
parSim.Opts.ca       = (parSim.OP.ca_ivc:parSim.Opts.deltaPhi:...
    parSim.OP.ca_evo);

%% Create acados Functions and Run Solver
% either solves for the optimization problem specified above ('SQP') or
% builds the s-func. necessary for the simulink implementation ('Build')

% if script is used to build s-functions for simulink implementation, only
% a single SQP-step is used, before re-applying the newly obtained solution
% if the optimization problem is to be solved with a matlab script, the 
% number of SQP-steps can be increased
switch parOpt.Solver
    case 'Build'
        parOpt.SQP.nSQPmax = 1;
    case 'SQP'
        parOpt.SQP.nSQPmax = 100;
end

% Extend optimization range until EVO for acados framework
parOpt.caAcados = parOpt.ca(end)+parOpt.deltaPhiAcados : ...
    parOpt.deltaPhiAcados : parSim.OP.ca_evo;
parOpt.nInt     = [repmat(parOpt.nInt,1,numel(parOpt.ca)-1), ...
    repmat(parOpt.nIntAcados,size(parOpt.caAcados))];
parOpt.ca       = [parOpt.ca, parOpt.ca(end)+parOpt.deltaPhiAcados : ...
    parOpt.deltaPhiAcados : parSim.OP.ca_evo];
switch parOpt.Solver
    case 'SQP' 
        ocp          = createAcadosFunctions_InjOpt(parOpt,parModel);
        uOpt         = runSQPacados_InjOpt(ocp,parModel,parSim,parOpt);
        
    case 'Build'
		[ocp,fcn]    = createAcadosFunctions_InjOpt(parOpt,parModel);
        runBuildAcados_InjOpt(ocp,fcn,parModel,parSim,parOpt);
		return
    otherwise
        error('No valid solver choice.');
end

%% Simulation With optimal Combustion
% Define Model Input and Initial Conditions
caSim  = parSim.Opts.ca;
x0Sim  = [parSim.OP.pInt,0,0]';
uSim   = repmat(uOpt,1,numel(caSim));

% Include NOx model
if parOpt.enNOx
    x0Sim = [x0Sim;parSim.OP.thetaIvc;0];
end

% Run simulation
simOpt = CompleteSimulation(caSim,x0Sim,uSim,parSim,parModel);

% Plot Results
fig2 = figure(2);
plotSimResultsVariation(fig2,simOpt,parSim);

fig3 = figure(3); clf
plotSimResultsDetails(fig3,simOpt,parSim,parOpt)

% run new optimization parameters with the following scrtip
% run testNewCase.m