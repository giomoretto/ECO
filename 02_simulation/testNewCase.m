% run this script, after the script 
%   'main_InjOpt.m', with the settings
%   parOpt.Solver = 'SQP'
% with the parameters below, you can change the bound values for the
% combustion, or the values of the operating point


% Define (new) reference values
parOpt.Reference.IMEP    = 8e5; % [Pa]
parOpt.Reference.pMax    = 150e5; % [Pa]
parOpt.Reference.dpMax   = 5e5; % [Pa/degCA]
parOpt.Reference.Tmin    = 0+273; % [K]
parOpt.Reference.cNOx    = 1e4; % [ppm]
parOpt.Reference.CoCmax  = 20; % [degCA]
parOpt.Reference.PhiMax  = 1/1.3; % [-]

% define (new) operating point
OP.Ne     = 2000/60; % [1/s];
OP.p_im   = 1.2e5; % [Pa] 
OP.T_im   = 273.15+45; % [K]
OP.x_bg   = 0.1; % [-]
OP.p_rail = 1000e5; % [Pa]

% calculate new operating point Parameters
parSim.OP = parOPDef(parModel,OP,'Num');

% Solve ocp, object ocp previousyl created in main_InjOpt.m
uOpt      = runSQPacados_InjOpt(ocp,parModel,parSim,parOpt);

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


