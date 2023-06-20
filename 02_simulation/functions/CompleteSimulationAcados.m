function simOpt = CompleteSimulationAcados(caSim,x0,u0,parSim,parModel)
%% Remove old folder (necessary if integrator is changed)
buildFolder       = what('06_build');
acadosBuildFolder = fullfile(buildFolder.path,'acadosBuild','Integrator');
if isfolder(acadosBuildFolder)
    rmpath(genpath(acadosBuildFolder));
    rmdir(acadosBuildFolder,'s');
    pause(0.5);
end

%% arguments
compile_interface = 'auto';
sim_method        = 'irk'; % erk, irk, irk_gnsf
num_stages        = 4;     % Runge-Kutta int. stages: (1) RK1, (4) RK4
num_steps         = 1;     % Intermediate RK4 steps without shooting node
newton_iter       = 5;     % for implicit integrators

% simulation parameters
nInputs           = parModel.nInputs;
nStates           = parModel.nStates;
N_sim             = numel(caSim);
h                 = caSim(2)-caSim(1);

%% acados sim model
import casadi.*
sim_model         = acados_sim_model();
sim_model.set('T', h);

%% define model dynamics
% Create symbolic model expressions with casADi
  % States (consisting of actual states and inputs, scaled)
    uxSym       = SX.sym('x',nInputs+nStates+1,1);
    uxdotSym    = SX.sym('xdot',nInputs+nStates+1,1);
  % Extract inputs and actual states
    u           = uxSym(1:nInputs);
    x           = uxSym(nInputs+1:nInputs+nStates);
    ca          = uxSym(nInputs+nStates+1);
  % System dynamics
    xdot      = CompleteModel(x,u,ca,parModel,parSim.OP,'Sym');
  % Get explicit and implicit derivative expressions for ux vector
    uxdot         = [zeros(nInputs,1); xdot; 1];
    uxdot_impl    = uxdot - uxdotSym;
  % Write to model object: symbolics
    sim_model.set('sym_x', uxSym);
    sim_model.set('sym_xdot', uxdotSym);
  % Write to model object: dynamics
    if (strcmp(sim_method, 'erk'))
        sim_model.set('dyn_type', 'explicit');
        sim_model.set('dyn_expr_f', uxdot);
    else % irk irk_gnsf
        sim_model.set('dyn_type', 'implicit');
        sim_model.set('dyn_expr_f', uxdot_impl);
    end

%% acados sim opts
sim_opts          = acados_sim_opts();

sim_opts.set('compile_interface', compile_interface);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('newton_iter', newton_iter); % for implicit intgrators
sim_opts.set('method', sim_method);
if (strcmp(sim_method, 'irk_gnsf'))
	sim_opts.set('gnsf_detect_struct', 'true');
end

% Set output direction
buildFolder       = what('06_build');
sim_opts.set('output_dir',fullfile(buildFolder.path,'acadosBuild','Integrator'));

%% create integrator
sim               = acados_sim(sim_model, sim_opts);

%% simulate system in loop
ux0               = [u0; x0; caSim(1)];
x_sim             = zeros(nInputs+nStates+1, N_sim);
xdot_sim          = zeros(nStates, N_sim);
y_sim             = zeros(parModel.nOutputs, N_sim);
x_sim(:,1)        = ux0;
[xdotT,yT]        = CompleteModel(x0,u0,caSim(1), ...
                        parModel,parSim.OP,'Num');
xdot_sim(:,1)     = xdotT;
y_sim(:,1)        = yT;
for ii=1:N_sim-1
	
	% set initial state
	sim.set('x', x_sim(:,ii));

    % initialize implicit integrator
    if (strcmp(sim_method, 'irk'))
        sim.set('xdot', zeros(nInputs+nStates+1,1));
    elseif (strcmp(sim_method, 'irk_gnsf'))
        n_out     = sim.model_struct.dim_gnsf_nout;
        sim.set('phi_guess', zeros(n_out,1));
    end

	% solve
	sim.solve();

	% get simulated state
	x_sim(:,ii+1)    = sim.get('xn');
    
	% get simulated state derivative and output
    [xdotT,yT]       = CompleteModel(x_sim(nInputs+1:nInputs+nStates,ii+1),u0,caSim(ii+1), ...
                        parModel,parSim.OP,'Num');
    xdot_sim(:,ii+1) = xdotT;
    y_sim(:,ii+1)    = yT;

end

% Create output struct
simOpt.ca         = caSim;
simOpt.x          = x_sim(nInputs+1:nInputs+nStates,:);
simOpt.u          = x_sim(1:nInputs,:);
simOpt.y          = y_sim;
simOpt.xdot       = xdot_sim;

end
